"""Rebuild Figures 1-3 from the supplied plotting tables; Figure 3 combines detection states and detected-subset abundance.
Writes editable SVG, 600-dpi RGB PNG and CMYK TIFF; never writes a PDF.
No statistical model is refitted and no research data are reselected.
"""
from pathlib import Path
import sys, os, argparse, hashlib, json, warnings, shutil
sys.stdout.reconfigure(encoding='utf-8')
HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'mpl_config'))
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
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,default=HERE/'rebuilt');ap.add_argument('--source-dir',type=Path,default=HERE/'source_data');ap.add_argument('--font-latin',type=Path);ap.add_argument('--font-cjk',type=Path);args=ap.parse_args()
DEST=args.output_dir.resolve();DEST.mkdir(parents=True,exist_ok=True)
DATA=args.source_dir.resolve()
VEC=DEST/'vector';PREV=DEST/'preview';OUT=DEST/'png';TIFF=DEST/'tiff'
for p in [VEC,PREV,OUT,TIFF]:p.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
T=['Lacticaseibacillus','Lactobacillus','Leuconostoc','Limosilactobacillus','Ligilactobacillus','Pediococcus']
CN=['乳酪杆菌属','乳杆菌属','明串珠菌属','黏液乳杆菌属','联合乳杆菌属','片球菌属']
BG=['Anaerococcus','Finegoldia','Peptoniphilus_A'];BCN=['厌氧球菌属','芬戈尔德菌属','嗜胨菌属']
AUDIT=[];TEXT_AUDIT={};EXPORTS={};PLOT_CHECKS=[]
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
C={'blue':'#0072B2','orange':'#D55E00','purple':'#7966A0','grey':'#92979D','grid':'#E9ECEF','flowblue':'#46545E','read1':'#46545E','read3':'#89939B'}
def figure(height):return plt.figure(figsize=(15/2.54,height/2.54),dpi=150)
def ftext(f,x,y,t,**kw):return f.text(x,y,t,fontsize=7,color='black',**kw)
def clean(ax,ylabels=None):
    for k in ['top','right','left']:ax.spines[k].set_visible(False)
    ax.tick_params(axis='y',length=0,pad=5)
    if ylabels is not None:ax.set_yticks(np.arange(len(ylabels)));ax.set_yticklabels(ylabels)
def err(ax,r,y,est='OR',low='low',high='high',color=C['blue'],marker='s',size=3.8,panel='',scale=1,filled=True):
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

def save(f,stem):
    f.canvas.draw();rend=f.canvas.get_renderer();labels=[];outside=[]
    for t in f.findobj(match=Text):
        if not t.get_visible() or not t.get_text():continue
        assert t.get_fontsize()==7,(stem,t.get_text(),t.get_fontsize())
        assert colors.to_rgba(t.get_color())==colors.to_rgba('black')
        bb=t.get_window_extent(rend);labels.append({'text':t.get_text(),'font_family':t.get_fontfamily(),'size_pt':t.get_fontsize(),'color':'black','bounds_pixels':list(bb.bounds)})
        if bb.x0<-.5 or bb.y0<-.5 or bb.x1>f.bbox.x1+.5 or bb.y1>f.bbox.y1+.5:outside.append(t.get_text())
    assert not outside,(stem,outside)
    TEXT_AUDIT[stem]={'labels':labels,'outside_canvas':outside}
    f.savefig(VEC/f'{stem}.svg')
    f.savefig(OUT/f'{stem}.png',dpi=600)
    f.savefig(PREV/f'{stem}.png',dpi=200)
    im=Image.open(OUT/f'{stem}.png').convert('RGB')
    im.save(OUT/f'{stem}.png',dpi=(600,600),icc_profile=SRGB.tobytes())
    assert im.width==3543,im.size
    cmyk=ImageCms.applyTransform(im,TRANSFORM)
    cmyk.save(TIFF/f'{stem}.tiff',compression='tiff_lzw',dpi=(600,600),icc_profile=CMYK.tobytes())
    png=Image.open(OUT/f'{stem}.png');tif=Image.open(TIFF/f'{stem}.tiff')
    assert abs(png.info['dpi'][0]-600)<.1 and tif.mode=='CMYK' and tif.info.get('icc_profile')==CMYK.tobytes()
    svg=(VEC/f'{stem}.svg').read_text(encoding='utf-8')
    assert '<text' in svg and '<image' not in svg,'SVG must retain live text and no embedded raster'
    proof=ImageCms.profileToProfile(cmyk,CMYK,SRGB,outputMode='RGB',renderingIntent=ImageCms.Intent.RELATIVE_COLORIMETRIC)
    proof.thumbnail((1600,1600));proof.save(PREV/f'{stem}_CMYK校样.png')
    EXPORTS[stem]={'pixels':im.size,'dpi':png.info['dpi'],'width_cm':15,'height_cm':float(f.get_figheight()*2.54),'fonts':names,'font_size_pt':7,'svg_editable_text':True,'svg_raster_images':0,'png_color':'RGB','tiff_color':tif.mode,'ICC':ICC_INFO,'files':{f'{sub}/{stem}.{ext}':sha(DEST/sub/f'{stem}.{ext}') for sub,ext in [('vector','svg'),('png','png'),('tiff','tiff')]}}
    im.close();png.close();tif.close();cmyk.close();proof.close();plt.close(f)

def fig1():
    f=figure(10);ax=f.add_axes([0,0,1,1]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    ftext(f,.024,.965,'(a) 研究对象筛选',va='top');ftext(f,.525,.965,'(b) 六属检出率',va='top')
    flow=json.loads((DATA/'fig1_flow.json').read_text(encoding='utf-8'))
    def box(x,y,w,h,text,blue=False):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.003,rounding_size=.003',edgecolor=C['flowblue'] if blue else '#B4BBC0',facecolor='#F1F3F4' if blue else 'white',lw=.7 if blue else .6))
        ax.text(x+w/2,y+h/2,text,ha='center',va='center',linespacing=1.15)
    def arrow(x0,y0,x1,y1,col='black'):ax.annotate('',xy=(x1,y1),xytext=(x0,y0),arrowprops=dict(arrowstyle='->',color=col,lw=.65,shrinkA=0,shrinkB=0,mutation_scale=7))
    x=.04;w=.225;h=.078;ys=[.80,.66,.52,.38,.24]
    texts=[f'粪便记录\n{flow["stool_records"]}份',f'每人首份粪便\n{flow["first_samples"]}人',f'年龄18～100岁\n{flow["adult_samples"]}人',f'关键变量完整\n{flow["complete_samples"]}人',f'基础样本描述\n{flow["descriptive_samples"]}人']
    exclusions=[f'其余粪便记录\n{flow["other_stools"]}份',f'年龄不符或缺失\n{flow["age_excluded"]}人',f'关键变量缺失\n{flow["missing_excluded"]}人',f'总读数<5000\n{flow["depth_excluded"]}人']
    cx=x+w/2
    for j,(y,t) in enumerate(zip(ys,texts)):
        box(x,y,w,h,t,blue=j==4)
        if j<4:
            arrow(cx,y,cx,ys[j+1]+h);ey=(y+ys[j+1]+h)/2
            arrow(cx,ey,.285,ey,'#89939B');ax.text(.302,ey,exclusions[j],ha='left',va='center',linespacing=1.12)
    box(x,.10,w,h,f'主关联分析\n{flow["association_samples"]}人',True)
    arrow(cx,ys[-1],cx,.178)
    arrow(cx,.209,.285,.209,'#89939B');ax.text(.302,.209,f'其他或未知性别\n{flow["sex_excluded"]}人',ha='left',va='center',linespacing=1.12)
    d=pd.read_csv(DATA/'fig1_detection.csv');p=f.add_axes([.649,.25,.205,.59]);p.set_xlim(0,21);p.set_ylim(5.5,-.5);clean(p,CN);p.set_xticks([0,5,10,15,20]);p.set_xlabel('检出率(%)',labelpad=5)
    for j,t in enumerate(T):
        p.hlines(j,0,21,color=C['grid'],lw=.55,zorder=0)
        for threshold,color,marker in [(1,C['read1'],'o'),(3,C['read3'],'s')]:
            r=d[(d.taxon==t)&(d.threshold==threshold)].iloc[0].to_dict();hnd=p.plot(r['percent'],j,marker=marker,color=color,markerfacecolor=color if threshold==1 else 'white',markeredgewidth=.85,ms=3.8,lw=0,zorder=3)[0]
            assert hnd.get_xdata()[0]==r['percent'];record('1b',r,['detected','n_total','percent','threshold'])
            p.text(1.20 if threshold==1 else 1.52,j,str(int(r['detected'])),transform=p.get_yaxis_transform(),ha='center',va='center',clip_on=False)
    p.text(1.36,1.13,'检出人数',transform=p.transAxes,ha='center',va='center')
    for xp,lab in [(1.20,'≥1条'),(1.52,'≥3条')]:p.text(xp,1.055,lab,transform=p.transAxes,ha='center',va='center')
    handles=[Line2D([],[],color=C['read1'],marker='o',linestyle='',markersize=3.8,label='≥1条读段'),Line2D([],[],color=C['read3'],marker='s',markerfacecolor='white',markeredgewidth=.85,linestyle='',markersize=3.8,label='≥3条读段')]
    f.legend(handles=handles,loc='lower center',bbox_to_anchor=(.79,.09),ncol=2,frameon=False,columnspacing=1.8,handletextpad=.7)
    save(f,'图1_分析人群与乳酸菌检出')

def fig2():
    f=figure(9.4);ftext(f,.012,.969,'(a) 六属背景关联比较',va='top');ftext(f,.49,.969,'(b) 基础调整后的检出关联',va='top');ftext(f,.49,.495,'(c) 饮食调整后的属间差异',va='top')
    a=f.add_axes([.151,.285,.283,.465]);a.set_xlim(-.5,4.5);a.set_ylim(4.5,-.5);a.set_aspect('equal');a.set_xticks(range(5));a.set_xticklabels(CN[:-1],rotation=48,ha='right',rotation_mode='anchor');a.set_yticks(range(5));a.set_yticklabels(CN[1:]);a.tick_params(length=0,pad=5)
    for sp in a.spines.values():sp.set_visible(False)
    cmap=colors.LinearSegmentedColormap.from_list('wald',['#F2F4F5','#A3AEB7']);norm=colors.Normalize(0,90)
    rows=table('C04')
    for i in range(1,6):
        for j in range(i):
            r=select(rows,taxon_A=T[j],taxon_B=T[i]);val=r['chi2'];sig=r['q_BH_15']<.05
            a.add_patch(Rectangle((j-.5,i-1-.5),1,1,facecolor=cmap(norm(val)),edgecolor='white',lw=.9))
            a.text(j,i-1,f'{val:.1f}'+('*' if sig else ''),ha='center',va='center',fontweight='bold' if (i,j)==(1,0) else 'normal');record('2a',r,['chi2','q_BH_15'])
    a.add_patch(Rectangle((-.5,-.5),1,1,facecolor='none',edgecolor='#41484D',lw=.9,zorder=4))
    cax=f.add_axes([.278,.837,.156,.016]);cb=f.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal');cb.set_ticks([0,30,60,90]);cb.solids.set_rasterized(False);cax.set_title(r'Wald $\chi^{2}$',pad=7);ftext(f,.294,.735,r'* $q$ < 0.05')
    b=f.add_axes([.635,.64,.345,.245]);clean(b,BCN);b.set_ylim(2.5,-.5);b.set_xlim(.75,1.8);b.set_xticks(np.arange(.8,1.81,.2));b.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'));b.set_xlabel('检出优势比(OR)',labelpad=5);b.axvline(1,color='#808080',ls=(0,(3,3)),lw=.7)
    rr=table('C13')
    for i,bg in enumerate(BG):
        b.axhline(i,color=C['grid'],lw=.45,zorder=0)
        for t,off,col,m in [(T[0],-.14,C['blue'],'s'),(T[1],.14,C['orange'],'o')]:
            r=select(rr,scenario='original_reads1',feature='microbe::Peptoniphilaceae|'+bg,taxon=t);err(b,r,i+off,low='CI_low',high='CI_high',color=col,marker=m,panel='2b')
    handles=[Line2D([],[],color=C['blue'],marker='s',lw=.8,markersize=3.6,label=CN[0]),Line2D([],[],color=C['orange'],marker='o',lw=.8,markersize=3.6,label=CN[1])]
    f.legend(handles=handles,loc='upper left',bbox_to_anchor=(.49,.935),ncol=2,frameon=False,columnspacing=1.1,handletextpad=.6,handlelength=1.1)
    c=f.add_axes([.635,.17,.345,.245]);clean(c,BCN);c.set_ylim(2.5,-.5);c.set_xlim(.5,1.04);c.set_xticks(np.arange(.5,1.01,.1));c.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'));c.axvline(1,color='#808080',ls=(0,(3,3)),lw=.7);c.set_xlabel('乳酪杆菌属 / 乳杆菌属 OR 比',labelpad=5)
    rr=table('P02')
    for i,bg in enumerate(BG):
        c.axhline(i,color=C['grid'],lw=.45,zorder=0)
        err(c,select(rr,spec='food10_nutrient10',background=bg),i,est='ratio_OR',color=C['purple'],marker='D',panel='2c')
    save(f,'图2_属间比较与饮食调整')

def fig3_combined():
    """Main Figure 3: four detection states above two detected-subset abundance panels.

    All models, rows, estimates, confidence intervals and axis scales are inherited.
    The two groups explicitly distinguish the full cohort from detected subsets.
    """
    f=figure(11)
    ftext(f,.024,.975,'(a) 四种检出状态（2748人）',va='top')
    hs=[Line2D([],[],color=C['grey'],marker='s',mfc='white',lw=.8,markersize=3.3,label='低背景（第10百分位）'),Line2D([],[],color=C['purple'],marker='o',lw=.8,markersize=3.3,label='高背景（第90百分位）')]
    f.legend(handles=hs,loc='upper center',bbox_to_anchor=(.66,.989),ncol=2,frameon=False,columnspacing=1.5,handlelength=1.1,handletextpad=.6)
    rr=table('S02');titles=['两属\n均未检出','乳酪杆菌属\n单独检出','乳杆菌属\n单独检出','两属\n共同检出']
    for j in range(4):
        a=f.add_axes([.177+j*.207,.604,.186,.224]);clean(a,BCN if j==0 else None);a.set_yticks(range(3));a.set_yticklabels(BCN if j==0 else []);a.set_ylim(2.5,-.5);a.set_xlim(0,80 if j==0 else 25);a.set_xticks([0,40,80] if j==0 else [0,10,20]);a.set_title(titles[j],pad=6,linespacing=1.15)
        for x in ([0,40,80] if j==0 else [0,10,20]):a.axvline(x,color=C['grid'],lw=.5,zorder=0)
        for i,bg in enumerate(BG):
            for lev,off,col,m,size in [('p10',-.14,C['grey'],'s',3.3),('p90',.14,C['purple'],'o',3.3)]:
                r=select(rr,threshold=1,background=bg,level=lev,state=j);err(a,r,i+off,est='probability',color=col,marker=m,size=size,panel='3a',scale=100,filled=lev=='p90')
    ftext(f,.58,.535,'调整后检出概率(%)',ha='center',va='top')

    ftext(f,.024,.475,'(b) 检出者相对丰度',va='top')
    rr=table('A01')
    for j,t in enumerate(T[:2]):
        a=f.add_axes([.177+j*.414,.117,.393,.227]);clean(a,BCN if j==0 else None);a.set_yticks(range(3));a.set_yticklabels(BCN if j==0 else []);a.set_ylim(2.5,-.5);a.set_xlim(.75,1.75);a.set_xticks(np.arange(.8,1.61,.2));a.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'));a.axvline(1,color='#808080',ls=(0,(3,3)),lw=.7)
        n=select(rr,threshold=1,background=BG[0],target=t)['n'];a.set_title(f'{CN[j]}（{n}人）',pad=10)
        for i,bg in enumerate(BG):
            a.axhline(i,color=C['grid'],lw=.5,zorder=0);r=select(rr,threshold=1,background=bg,target=t);err(a,r,i,est='ratio_geomean',color=C['blue'] if j==0 else C['orange'],marker='s' if j==0 else 'o',panel='3b')
    ftext(f,.58,.045,'相对丰度几何均值比',ha='center',va='bottom')
    save(f,'图3_背景菌与检出状态及相对丰度')

with warnings.catch_warnings(record=True) as warning_list:
    warnings.simplefilter('always');fig1();fig2();fig3_combined()
    warns=sorted(set(str(w.message) for w in warning_list))
    assert not warns,warns
assert len(AUDIT)==66,len(AUDIT)
assert len(PLOT_CHECKS)==39,len(PLOT_CHECKS)
original=json.loads((HERE/'source_data/original_plotted_values.json').read_text(encoding='utf-8'))
# Both current Figure 3 groups retain their original audit panel identifiers.
legacy_panel={}
assert [(legacy_panel.get(x['panel'],x['panel']),x['identity'],x['values']) for x in AUDIT]==[(x['panel'],x['identity'],x['values']) for x in original]
(DEST/'numeric_plot_audit.json').write_text(json.dumps({'all_plotted_rows':AUDIT,'artist_ci_checks':PLOT_CHECKS,'n_plotted_rows':len(AUDIT),'all_artist_checks_pass':True,'exact_original_plot_value_equality':True,'statistical_models_refitted':False},ensure_ascii=False,indent=2),encoding='utf-8')
(DEST/'typography_audit.json').write_text(json.dumps(TEXT_AUDIT,ensure_ascii=False,indent=2),encoding='utf-8')
(DEST/'export_audit.json').write_text(json.dumps({'exports':EXPORTS,'warnings':warns,'matplotlib_version':matplotlib.__version__,'PDF_files_generated':0},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'output':str(DEST),'figures':len(EXPORTS),'plotted_rows':len(AUDIT),'CI_artist_checks':len(PLOT_CHECKS),'exact_original_value_equality':True,'PDF_files_generated':0},ensure_ascii=False,indent=2))
