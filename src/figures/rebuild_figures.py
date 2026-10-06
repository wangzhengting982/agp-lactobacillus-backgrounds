"""Render Figures 1-3 from explicitly prepared, provenance-tracked source tables.
Run via scripts/run_figures.py after statistical runs have completed.
Requires Python 3, numpy, pandas, matplotlib, pypdf and Pillow with LittleCMS2.
Uses installed Windows Arial and YouYuan fonts; profiles are bundled under color_profiles.
No model fitting, data reselection, numeric changes or edits to submission documents.
"""
from pathlib import Path
import sys, os, argparse, hashlib, json, warnings, shutil
sys.stdout.reconfigure(encoding='utf-8')
HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE.parents[1]/'runs/figures/mpl_config'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, colors, ticker
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image, ImageCms
from pypdf import PdfReader, PdfWriter
from pypdf.generic import ContentStream, NameObject, FloatObject, NumberObject, DictionaryObject, ArrayObject, DecodedStreamObject, TextStringObject
ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--source-dir',type=Path,required=True);ap.add_argument('--font-latin',type=Path);ap.add_argument('--font-cjk',type=Path);args=ap.parse_args()
DEST=args.output_dir.resolve();DEST.mkdir(parents=True,exist_ok=True)
DATA=args.source_dir.resolve()
VEC=DEST/'vector';PREV=DEST/'preview';OUT=DEST/'png';TIFF=DEST/'tiff'
RGBPDF=DEST/'source_rgb_pdf'
for p in [VEC,PREV,OUT,TIFF,RGBPDF]:p.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
T=['Lacticaseibacillus','Lactobacillus','Leuconostoc','Limosilactobacillus','Ligilactobacillus','Pediococcus']
CN=['乳酪杆菌属','乳杆菌属','明串珠菌属','黏液乳杆菌属','联合乳杆菌属','片球菌属']
BG=['Anaerococcus','Finegoldia','Peptoniphilus_A'];BCN=['厌氧球菌属','芬戈尔德菌属','嗜胨菌属']
AUDIT=[];TEXT_AUDIT={};EXPORTS={};PLOT_CHECKS=[];PDF_COLOR_AUDIT={}
TABLES=json.loads((DATA/'plot_tables.json').read_text(encoding='utf-8'))
def table(key):return TABLES[key]
def select(rows,**filters):
    rr=[r for r in rows if all(r[k]==v for k,v in filters.items())]
    assert len(rr)==1,(filters,len(rr));return rr[0]
def record(panel,row,fields,plotted=None):
    d={'panel':panel,'source':row.get('_source','verified cohort inputs'),'source_row':row.get('_row'),'identity':{k:row[k] for k in ['taxon','taxon_A','taxon_B','target','background','feature','exposure','state','level','threshold','spec','scenario'] if k in row},'values':{k:row[k] for k in fields}}
    if plotted is not None:d['plotted']=plotted
    AUDIT.append(d)

windows_fonts=Path(os.environ.get('WINDIR','C:/Windows'))/'Fonts'
fontpaths=[str(args.font_latin or windows_fonts/'arial.ttf'),str(args.font_cjk or windows_fonts/'SIMYOU.TTF')]
for p in fontpaths: assert Path(p).exists();font_manager.fontManager.addfont(p)
names=[font_manager.FontProperties(fname=p).get_name() for p in fontpaths]
plt.rcParams.update({'font.family':names,'font.size':7,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.labelsize':7,'axes.titlesize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':0.65,'xtick.major.width':0.65,'xtick.major.size':2.5,'ytick.major.width':0.65,'ytick.major.size':2.5,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False,'savefig.facecolor':'white','mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold','mathtext.fallback':None})
C={'blue':'#3387A8','orange':'#DC813D','purple':'#7966A0','grey':'#92979D','grid':'#E6E9EB','flowblue':'#007CB5','read1':'#007BB5','read3':'#D85D00'}
def figure(height):return plt.figure(figsize=(15/2.54,height/2.54),dpi=150)
def ftext(f,x,y,t,**kw):return f.text(x,y,t,fontsize=7,color='black',**kw)
def clean(ax,ylabels=None):
    for k in ['top','right','left']:ax.spines[k].set_visible(False)
    ax.tick_params(axis='y',length=0,pad=5)
    if ylabels is not None:ax.set_yticks(np.arange(len(ylabels)));ax.set_yticklabels(ylabels)
def err(ax,r,y,est='OR',low='low',high='high',color=C['blue'],marker='s',size=3.6,panel='',scale=1,filled=True):
    x,lo,hi=[float(r[k])*scale for k in [est,low,high]];assert lo<=x<=hi
    h=ax.errorbar(x,y,xerr=[[x-lo],[hi-x]],fmt=marker,markersize=size,mfc=color if filled else 'white',mec=color,mew=.7,color=color,elinewidth=.8,capsize=2.0,capthick=.7,zorder=3)
    xv=float(h.lines[0].get_xdata()[0]);seg=h.lines[2][0].get_segments()[0]
    assert np.isclose(xv,x,rtol=0,atol=1e-13) and np.allclose(seg[:,0],[lo,hi],rtol=0,atol=1e-13)
    PLOT_CHECKS.append({'panel':panel,'estimate':x,'low':lo,'high':hi,'artist_equal_source':True})
    record(panel,r,[est,low,high],{'estimate':x,'low':lo,'high':hi,'scale':scale});return h

SRGB_FILE=HERE/'color_profiles/sRGB_IEC61966-2-1.icm'
CMYK_FILE=HERE/'color_profiles/Agfa_SWOP_Standard.icm'
SRGB=ImageCms.getOpenProfile(str(SRGB_FILE));CMYK=ImageCms.getOpenProfile(str(CMYK_FILE))
TRANSFORM=ImageCms.buildTransformFromOpenProfiles(SRGB,CMYK,'RGB','CMYK',renderingIntent=ImageCms.Intent.RELATIVE_COLORIMETRIC,flags=ImageCms.Flags.BLACKPOINTCOMPENSATION)
ICC_INFO={'source_profile':ImageCms.getProfileDescription(SRGB).strip(),'destination_profile':ImageCms.getProfileDescription(CMYK).strip(),'source_space':SRGB.profile.xcolor_space.strip(),'destination_space':CMYK.profile.xcolor_space.strip(),'intent':'Relative colorimetric','black_point_compensation':True,'source_sha256':sha(SRGB_FILE),'destination_sha256':sha(CMYK_FILE),'conversion_engine':'Pillow ImageCms / LittleCMS2','profile_origin':'Windows System32/spool/drivers/color'}

def pdf_to_cmyk(source,destination):
    """Convert only verified vector color operators; preserve text and geometry."""
    reader=PdfReader(source);page=reader.pages[0];resources=page['/Resources']
    assert len(reader.pages)==1
    assert all(k not in resources or len(resources[k].get_object())==0 for k in ['/XObject','/Shading','/Pattern'])
    contents=ContentStream(page.get_contents(),reader)
    allowed={b'rg',b'RG',b'g',b'G',b'cs',b'CS'}
    assert not any(op in {b'sc',b'SC',b'scn',b'SCN',b'sh',b'Do',b'BI',b'ID',b'EI'} for _,op in contents.operations)
    before_ops=[(args,op) for args,op in contents.operations if op not in allowed]
    initial_text=page.extract_text();before_fonts={k:v.get_object()['/DescendantFonts'][0].get_object()['/FontDescriptor'].get_object()['/FontFile2'].get_object().get_data() for k,v in resources['/Font'].items()}
    rgb_cache={};changed={'RGB_fill':0,'RGB_stroke':0,'Gray_fill':0,'Gray_stroke':0,'space':0}
    new=[]
    for args,op in contents.operations:
        if op in {b'rg',b'RG'}:
            rgb=tuple(max(0,min(255,round(float(x)*255))) for x in args)
            if rgb not in rgb_cache:rgb_cache[rgb]=ImageCms.applyTransform(Image.new('RGB',(1,1),rgb),TRANSFORM).getpixel((0,0))
            new.append(([FloatObject(x/255) for x in rgb_cache[rgb]],b'k' if op==b'rg' else b'K'));changed['RGB_fill' if op==b'rg' else 'RGB_stroke']+=1
        elif op in {b'g',b'G'}:
            # DeviceGray represents neutral coverage; retain it as K-only CMYK.
            new.append(([FloatObject(0),FloatObject(0),FloatObject(0),FloatObject(1-float(args[0]))],b'k' if op==b'g' else b'K'));changed['Gray_fill' if op==b'g' else 'Gray_stroke']+=1
        elif op in {b'cs',b'CS'}:
            assert str(args[0])=='/DeviceRGB',args
            new.append(([NameObject('/DeviceCMYK')],op));changed['space']+=1
        else:new.append((args,op))
    contents.operations=new
    assert before_ops==[(args,op) for args,op in new if op not in {b'k',b'K',b'cs',b'CS'}]
    writer=PdfWriter();writer.add_page(page);wp=writer.pages[0]
    wp[NameObject('/Contents')]=writer._add_object(contents)
    profile=DecodedStreamObject();profile.set_data(CMYK.tobytes());profile.update({NameObject('/N'):NumberObject(4),NameObject('/Alternate'):NameObject('/DeviceCMYK')});pref=writer._add_object(profile)
    wp['/Resources'][NameObject('/ColorSpace')]=DictionaryObject({NameObject('/DefaultCMYK'):ArrayObject([NameObject('/ICCBased'),pref])})
    intent=DictionaryObject({NameObject('/Type'):NameObject('/OutputIntent'),NameObject('/S'):NameObject('/GTS_PDFX'),NameObject('/OutputConditionIdentifier'):TextStringObject('Agfa SWOP Standard'),NameObject('/Info'):TextStringObject('Agfa SWOP Standard; ICC color conversion; no PDF/X conformance claim'),NameObject('/DestOutputProfile'):pref})
    writer._root_object[NameObject('/OutputIntents')]=ArrayObject([writer._add_object(intent)])
    writer.add_metadata({'/Title':source.stem,'/Creator':'Matplotlib; pypdf vector CMYK conversion with LittleCMS2'})
    with destination.open('wb') as f:writer.write(f)
    final=PdfReader(destination);fp=final.pages[0];operations=ContentStream(fp.get_contents(),final).operations
    assert fp.extract_text()==initial_text
    for k,v in fp['/Resources']['/Font'].items():assert v.get_object()['/DescendantFonts'][0].get_object()['/FontDescriptor'].get_object()['/FontFile2'].get_object().get_data()==before_fonts[k]
    assert not any(op in {b'rg',b'RG',b'g',b'G',b'sc',b'SC',b'scn',b'SCN'} for _,op in operations)
    assert not any(op in {b'cs',b'CS'} and str(args[0])!='/DeviceCMYK' for args,op in operations)
    PDF_COLOR_AUDIT[source.stem]={'converted_operators':changed,'RGB_colors_converted':len(rgb_cache),'remaining_RGB_operators':0,'all_color_spaces_CMYK':True,'embedded_ICC_space':CMYK.profile.xcolor_space.strip(),'embedded_ICC_sha256':hashlib.sha256(CMYK.tobytes()).hexdigest(),'K_only_original_gray':True,'text_extraction_identical':True,'font_streams_identical':True,'geometry_operations_identical':True,'images':0,'shadings':0,'patterns':0,'PDF_X_conformance_claim':False}

def save(f,stem):
    f.canvas.draw();rend=f.canvas.get_renderer();labels=[];outside=[]
    for t in f.findobj(match=Text):
        if not t.get_visible() or not t.get_text():continue
        assert t.get_fontsize()==7,(stem,t.get_text(),t.get_fontsize())
        assert colors.to_rgba(t.get_color())==colors.to_rgba('black')
        bb=t.get_window_extent(rend);labels.append({'text':t.get_text(),'font_family':t.get_fontfamily(),'size_pt':t.get_fontsize(),'color':'black','bounds_pixels':list(bb.bounds)})
        if bb.x0<-.5 or bb.y0<-.5 or bb.x1>f.bbox.x1+.5 or bb.y1>f.bbox.y1+.5:outside.append(t.get_text())
    overlaps=[]
    # Disjoint text rectangles are guaranteed not to overlap. Rotated labels are
    # assessed visually, because their axis-aligned rectangles overestimate ink.
    for i,a in enumerate(labels):
        x,y,w,h=a['bounds_pixels']
        for b in labels[i+1:]:
            xx,yy,ww,hh=b['bounds_pixels']
            ix=min(x+w,xx+ww)-max(x,xx);iy=min(y+h,yy+hh)-max(y,yy)
            if ix>1 and iy>1 and ix*iy>0.08*min(w*h,ww*hh):overlaps.append([a['text'],b['text']])
    TEXT_AUDIT[stem]={'labels':labels,'outside_canvas':outside,'overlap_candidates':overlaps};assert not outside,(stem,outside)
    f.savefig(VEC/f'{stem}.svg')
    f.savefig(RGBPDF/f'{stem}.pdf')
    pdf_to_cmyk(RGBPDF/f'{stem}.pdf',VEC/f'{stem}.pdf')
    f.savefig(OUT/f'{stem}.png',dpi=600)
    f.savefig(PREV/f'{stem}.png',dpi=200)
    im=Image.open(OUT/f'{stem}.png').convert('RGB')
    im.save(OUT/f'{stem}.png',dpi=(600,600),icc_profile=SRGB.tobytes())
    assert im.width==3543,im.size
    # Explicit ICC conversion, not Image.convert('CMYK').
    cmyk=ImageCms.applyTransform(im,TRANSFORM)
    cmyk.save(TIFF/f'{stem}.tiff',compression='tiff_lzw',dpi=(600,600),icc_profile=CMYK.tobytes())
    png=Image.open(OUT/f'{stem}.png');tif=Image.open(TIFF/f'{stem}.tiff')
    assert abs(png.info['dpi'][0]-600)<.1 and tif.mode=='CMYK' and tif.info.get('icc_profile')==CMYK.tobytes()
    svg=(VEC/f'{stem}.svg').read_text(encoding='utf-8')
    assert '<text' in svg and '<image' not in svg,'SVG must be entirely vector with live text'
    # Roundtrip preview supports visual inspection of the CMYK interpretation.
    proof=ImageCms.profileToProfile(cmyk,CMYK,SRGB,outputMode='RGB',renderingIntent=ImageCms.Intent.RELATIVE_COLORIMETRIC)
    proof.thumbnail((1600,1600));proof.save(PREV/f'{stem}_CMYK校样.png')
    EXPORTS[stem]={'pixels':im.size,'dpi':png.info['dpi'],'vector_width_cm':15,'raster_width_cm':im.width/600*2.54,'height_cm':float(f.get_figheight()*2.54),'fonts':names,'font_size_pt':7,'font_color':'black','svg_editable_text':True,'svg_raster_images':0,'pdf_color':'CMYK vector with embedded ICC','png_color':'RGB','tiff_color':tif.mode,'ICC':ICC_INFO,'files':{f'{sub}/{stem}.{ext}':sha(DEST/sub/f'{stem}.{ext}') for sub,ext in [('vector','svg'),('vector','pdf'),('png','png'),('tiff','tiff')]}}
    im.close();png.close();tif.close();cmyk.close();proof.close();plt.close(f)
def fig1():
    f=figure(10);ax=f.add_axes([0,0,1,1]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    ftext(f,.024,.965,'(a) 研究对象筛选',va='top');ftext(f,.525,.965,'(b) 六属检出率',va='top')
    flow=json.loads((DATA/'fig1_flow.json').read_text(encoding='utf-8'))
    def box(x,y,w,h,text,blue=False):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.004,rounding_size=.004',edgecolor=C['flowblue'] if blue else '#A0A5A8',facecolor='#F1F5F7' if blue else 'white',lw=.65))
        ax.text(x+w/2,y+h/2,text,ha='center',va='center',linespacing=1.15)
    def arrow(x0,y0,x1,y1,col='black'):ax.annotate('',xy=(x1,y1),xytext=(x0,y0),arrowprops=dict(arrowstyle='->',color=col,lw=.65,shrinkA=0,shrinkB=0,mutation_scale=7))
    x=.035;w=.205;h=.078;ys=[.78,.662,.544,.426,.308]
    texts=[f'粪便记录\n{flow["stool_records"]}份',f'每人首份粪便\n{flow["first_samples"]}人',f'年龄18～100岁\n{flow["adult_samples"]}人',f'关键变量完整\n{flow["complete_samples"]}人',f'总读数≥5000\n{flow["descriptive_samples"]}人']
    exclusions=[f'其余粪便记录\n{flow["other_stools"]}份',f'年龄不符或缺失\n{flow["age_excluded"]}人',f'关键变量缺失\n{flow["missing_excluded"]}人',f'总读数<5000\n{flow["depth_excluded"]}人']
    cx=x+w/2
    for j,(y,t) in enumerate(zip(ys,texts)):
        box(x,y,w,h,t,blue=j==4)
        if j<4:
            arrow(cx,y,cx,ys[j+1]+h);ey=(y+ys[j+1]+h)/2
            arrow(cx,ey,.275,ey,'#808080');ax.text(.287,ey,exclusions[j],ha='left',va='center',linespacing=1.12)
    box(x,.098,w,h,f'主关联分析\n{flow["association_samples"]}人',True)
    box(.29,.098,.205,h,f'基础样本描述\n{flow["descriptive_samples"]}人')
    arrow(cx,ys[-1],cx,.176);ax.plot([cx,.3925,.3925],[.267,.267,.177],color='black',lw=.65)
    arrow(.3925,.197,.3925,.176)
    arrow(cx,.216,.248,.216,'#808080');ax.text(.259,.216,f'其他或未知性别\n{flow["sex_excluded"]}人',ha='left',va='center',linespacing=1.12)
    d=pd.read_csv(DATA/'fig1_detection.csv');p=f.add_axes([.649,.264,.205,.551]);p.set_xlim(0,21);p.set_ylim(5.5,-.5);clean(p,CN);p.set_xticks([0,5,10,15,20]);p.set_xlabel('检出率(%)',labelpad=5)
    for j,t in enumerate(T):
        p.hlines(j,0,21,color=C['grid'],lw=.55,zorder=0)
        for threshold,color,marker in [(1,C['read1'],'o'),(3,C['read3'],'s')]:
            r=d[(d.taxon==t)&(d.threshold==threshold)].iloc[0].to_dict();hnd=p.plot(r['percent'],j,marker=marker,color=color,ms=4,lw=0)[0]
            assert hnd.get_xdata()[0]==r['percent'];record('1b',r,['detected','n_total','percent','threshold'])
            p.text(1.20 if threshold==1 else 1.52,j,str(int(r['detected'])),transform=p.get_yaxis_transform(),ha='center',va='center',clip_on=False)
    p.text(1.36,1.22,'检出人数',transform=p.transAxes,ha='center',va='center')
    for xp,lab in [(1.20,'≥1条'),(1.52,'≥3条')]:p.text(xp,1.145,lab,transform=p.transAxes,ha='center',va='center')
    handles=[Line2D([],[],color=C['read1'],marker='o',linestyle='',markersize=4,label='≥1条读段'),Line2D([],[],color=C['read3'],marker='s',linestyle='',markersize=4,label='≥3条读段')]
    f.legend(handles=handles,loc='lower center',bbox_to_anchor=(.79,.065),ncol=2,frameon=False,columnspacing=1.8,handletextpad=.7)
    save(f,'图1_分析人群与乳酸菌检出')

def fig2():
    f=figure(10.8);ftext(f,.012,.969,'(a) 六属两两比较',va='top');ftext(f,.565,.969,'(b) 背景菌与两属检出的关联',va='top');ftext(f,.565,.479,'(c) 饮食调整后的属间差异',va='top')
    a=f.add_axes([.165,.29,.300,.465]);a.set_xlim(-.5,4.5);a.set_ylim(4.5,-.5);a.set_aspect('equal');a.set_xticks(range(5));a.set_xticklabels(CN[:-1],rotation=48,ha='right',rotation_mode='anchor');a.set_yticks(range(5));a.set_yticklabels(CN[1:]);a.tick_params(length=0,pad=5)
    for sp in a.spines.values():sp.set_visible(False)
    cmap=colors.LinearSegmentedColormap.from_list('wald',['#F0F4F6','#8CB7CA']);norm=colors.Normalize(0,90)
    rows=table('C04')
    for i in range(1,6):
        for j in range(i):
            r=select(rows,taxon_A=T[j],taxon_B=T[i]);val=r['chi2'];sig=r['q_BH_15']<.05
            a.add_patch(Rectangle((j-.5,i-1-.5),1,1,facecolor=cmap(norm(val)),edgecolor='white',lw=.9))
            a.text(j,i-1,f'{val:.1f}'+('*' if sig else ''),ha='center',va='center');record('2a',r,['chi2','q_BH_15'])
    cax=f.add_axes([.33,.83,.16,.016]);cb=f.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal');cb.set_ticks([0,30,60,90]);cb.solids.set_rasterized(False);cax.set_title(r'Wald $\chi^{2}$',pad=7);ftext(f,.321,.696,r'* $q$ < 0.05')
    b=f.add_axes([.72,.64,.255,.24]);clean(b,BCN);b.set_ylim(2.5,-.5);b.set_xlim(.75,1.8);b.set_xticks(np.arange(.8,1.81,.2));b.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'));b.set_xlabel('检出优势比(OR)',labelpad=5);b.axvline(1,color='#808080',ls=(0,(3,3)),lw=.7)
    rr=table('C13')
    for i,bg in enumerate(BG):
        for t,off,col,m in [(T[0],-.14,C['blue'],'s'),(T[1],.14,C['orange'],'o')]:
            r=select(rr,scenario='original_reads1',feature='microbe::Peptoniphilaceae|'+bg,taxon=t);err(b,r,i+off,low='CI_low',high='CI_high',color=col,marker=m,panel='2b')
    handles=[Line2D([],[],color=C['blue'],marker='s',lw=.8,markersize=3.6,label=CN[0]),Line2D([],[],color=C['orange'],marker='o',lw=.8,markersize=3.6,label=CN[1])]
    f.legend(handles=handles,loc='upper left',bbox_to_anchor=(.565,.935),ncol=2,frameon=False,columnspacing=1.1,handletextpad=.6,handlelength=1.1)
    c=f.add_axes([.72,.172,.255,.235]);clean(c,BCN);c.set_ylim(2.5,-.5);c.set_xlim(.5,1.04);c.set_xticks(np.arange(.5,1.01,.1));c.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'));c.axvline(1,color='#808080',ls=(0,(3,3)),lw=.7);c.set_xlabel('乳酪杆菌属 / 乳杆菌属 OR 比',labelpad=5)
    rr=table('P02')
    for i,bg in enumerate(BG):err(c,select(rr,spec='food10_nutrient10',background=bg),i,est='ratio_OR',color=C['purple'],marker='D',panel='2c')
    save(f,'图2_属间比较与饮食调整')

def fig3():
    f=figure(10);ftext(f,.012,.97,'(a) 两属的检出概率',va='top');ftext(f,.012,.465,'(b) 背景菌与目标属相对丰度的关联',va='top')
    hs=[Line2D([],[],color=C['grey'],marker='s',lw=.8,markersize=2.4,label='背景菌低水平'),Line2D([],[],color=C['purple'],marker='o',lw=.8,markersize=3.4,label='背景菌高水平')]
    f.legend(handles=hs,loc='upper right',bbox_to_anchor=(.982,.984),ncol=2,frameon=False,columnspacing=1.2,handlelength=1.1,handletextpad=.6)
    rr=table('S02');titles=['两属\n均未检出','乳酪杆菌属\n单独检出','乳杆菌属\n单独检出','两属\n共同检出']
    for j in range(4):
        a=f.add_axes([.177+j*.207,.605,.186,.25]);clean(a,BCN if j==0 else None);a.set_yticks(range(3));a.set_yticklabels(BCN if j==0 else []);a.set_ylim(2.5,-.5);a.set_xlim(0,80 if j==0 else 25);a.set_xticks([0,40,80] if j==0 else [0,10,20]);a.set_title(titles[j],pad=7,linespacing=1.1)
        for x in ([0,40,80] if j==0 else [0,10,20]):a.axvline(x,color=C['grid'],lw=.5,zorder=0)
        for i,bg in enumerate(BG):
            for lev,off,col,m,size in [('p10',-.14,C['grey'],'s',2.4),('p90',.14,C['purple'],'o',3.4)]:
                r=select(rr,threshold=1,background=bg,level=lev,state=j);err(a,r,i+off,est='probability',color=col,marker=m,size=size,panel='3a',scale=100)
    ftext(f,.58,.545,'调整后检出概率(%)',ha='center',va='top')
    rr=table('A01')
    for j,t in enumerate(T[:2]):
        a=f.add_axes([.177+j*.414,.105,.393,.25]);clean(a,BCN if j==0 else None);a.set_yticks(range(3));a.set_yticklabels(BCN if j==0 else []);a.set_ylim(2.5,-.5);a.set_xlim(.75,1.75);a.set_xticks(np.arange(.8,1.61,.2));a.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'));a.axvline(1,color='#808080',ls=(0,(3,3)),lw=.7)
        n=select(rr,threshold=1,background=BG[0],target=t)['n'];a.set_title(f'{CN[j]}({n}人)',pad=10)
        for i,bg in enumerate(BG):
            a.axhline(i,color=C['grid'],lw=.5,zorder=0);r=select(rr,threshold=1,background=bg,target=t);err(a,r,i,est='ratio_geomean',color=C['blue'] if j==0 else C['orange'],marker='s' if j==0 else 'o',panel='3b')
    ftext(f,.58,.024,'相对丰度几何均值比',ha='center',va='bottom');save(f,'图3_背景水平与检出及丰度')


with warnings.catch_warnings(record=True) as warning_list:
    warnings.simplefilter('always');fig1();fig2();fig3()
    warns=sorted(set(str(w.message) for w in warning_list))
    assert not warns,warns
assert len(AUDIT)==66,len(AUDIT)
assert len(PLOT_CHECKS)==39,len(PLOT_CHECKS)
original=json.loads((HERE/'source_data/original_plotted_values.json').read_text(encoding='utf-8'))
# Identity, plotted estimates, CI endpoints, percentages and heatmap statistics
# must exactly match the previous exported artwork. Only style and size change.
assert len(AUDIT)==len(original)
numeric_differences=[]
for new,old in zip(AUDIT,original):
    assert (new['panel'],new['identity'])==(old['panel'],old['identity'])
    assert new['values'].keys()==old['values'].keys()
    for key,value in new['values'].items():
        reference=old['values'][key]
        if isinstance(value,(float,int)) and isinstance(reference,(float,int)):
            assert np.isclose(value,reference,rtol=1e-7,atol=1e-9),(new['panel'],key,value,reference)
            if value!=reference:numeric_differences.append({'panel':new['panel'],'identity':new['identity'],'field':key,'rerun':value,'submitted':reference,'absolute_difference':abs(value-reference)})
        else:assert value==reference
(DEST/'numeric_plot_audit.json').write_text(json.dumps({'all_plotted_rows':AUDIT,'artist_ci_checks':PLOT_CHECKS,'n_plotted_rows':len(AUDIT),'all_artist_checks_pass':True,'exact_original_plot_value_equality':len(numeric_differences)==0,'all_original_plot_values_within_tolerance':True,'tolerance':{'relative':1e-7,'absolute':1e-9},'numeric_differences':numeric_differences},ensure_ascii=False,indent=2),encoding='utf-8')
(DEST/'typography_audit.json').write_text(json.dumps(TEXT_AUDIT,ensure_ascii=False,indent=2),encoding='utf-8')
(DEST/'pdf_color_audit.json').write_text(json.dumps(PDF_COLOR_AUDIT,ensure_ascii=False,indent=2),encoding='utf-8')
(DEST/'export_audit.json').write_text(json.dumps({'exports':EXPORTS,'warnings':warns,'matplotlib_version':matplotlib.__version__},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'output':str(DEST),'figures':len(EXPORTS),'plotted_rows':len(AUDIT),'CI_artist_checks':len(PLOT_CHECKS),'all_values_within_declared_tolerance':True,'exact_original_value_equality':len(numeric_differences)==0,'ICC':ICC_INFO,'overlap_candidates':{k:v['overlap_candidates'] for k,v in TEXT_AUDIT.items()}},ensure_ascii=False,indent=2))
