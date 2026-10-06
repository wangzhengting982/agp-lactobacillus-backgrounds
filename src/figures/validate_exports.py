"""Validate emitted vector and raster artwork; render final CMYK PDFs for review."""
from pathlib import Path
import json,io,sys,hashlib,argparse
from pypdf import PdfReader
from pypdf.generic import ContentStream
from PIL import Image
from lxml import etree
from fontTools.ttLib import TTFont
import pypdfium2 as pdfium
sys.stdout.reconfigure(encoding='utf-8')
ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
r=args.output_dir.resolve();audit={}
for pdf in sorted((r/'vector').glob('*.pdf')):
 reader=PdfReader(pdf);page=reader.pages[0];fonts=[]
 for k,v in page['/Resources']['/Font'].items():
  x=v.get_object();target=x['/DescendantFonts'][0].get_object();fd=target['/FontDescriptor'].get_object();fontbytes=fd['/FontFile2'].get_object().get_data();tt=TTFont(io.BytesIO(fontbytes))
  fonts.append({'resource':str(k),'basefont':str(x['/BaseFont']),'subtype':str(x['/Subtype']),'embedded_FontFile2':True,'FontName':str(fd['/FontName']),'embedded_post_italic_angle':float(tt['post'].italicAngle),'embedded_italic_flag':bool(tt['head'].macStyle&2),'embedded_font_sha256':hashlib.sha256(fontbytes).hexdigest()})
  assert any(n in str(x['/BaseFont']) for n in ['YouYuan','ArialMT','Arial-ItalicMT'])
 tree=etree.parse(str(pdf.with_suffix('.svg')));ns={'s':'http://www.w3.org/2000/svg'}
 math_nodes=[{'text':''.join(x.itertext()).strip(),'style':x.get('style')} for x in tree.xpath('//s:tspan',namespaces=ns)]
 assert not tree.xpath('//s:image',namespaces=ns)
 assert len(tree.xpath('//s:text',namespaces=ns))>0
 if pdf.stem.startswith('图2'):
  for char in ['q','χ']:assert any(x['text']==char and 'italic' in x['style'] and "'Arial'" in x['style'] for x in math_nodes)
  assert any(x['embedded_italic_flag'] and x['embedded_post_italic_angle']<0 for x in fonts)
 visible_text=''.join(tree.getroot().itertext());assert '（' not in visible_text and '）' not in visible_text
 png=Image.open(r/'png'/(pdf.stem+'.png'));tif=Image.open(r/'tiff'/(pdf.stem+'.tiff'))
 assert png.mode=='RGB' and tif.mode=='CMYK' and png.size==tif.size
 assert abs(float(page.mediabox.width)/72*2.54-15)<1e-6
 assert all(abs(float(d)-600)<.01 for d in png.info['dpi']) and all(float(d)==600 for d in tif.info['dpi'])
 ops=ContentStream(page.get_contents(),reader).operations
 assert not any(op in [b'rg',b'RG',b'g',b'G',b'sc',b'SC',b'scn',b'SCN'] for _,op in ops)
 assert page['/Resources']['/ColorSpace']['/DefaultCMYK'][0]=='/ICCBased'
 assert page['/Resources']['/ColorSpace']['/DefaultCMYK'][1].get_object()['/N']==4
 assert not any(op in [b'Do',b'sh'] for _,op in ops)
 doc=pdfium.PdfDocument(str(pdf));pp=doc[0];bitmap=pp.render(scale=200/72);im=bitmap.to_pil();im.save(r/'preview'/(pdf.stem+'_CMYK_PDF.png'));pp.close();doc.close()
 data={'page_cm':[float(page.mediabox.width)/72*2.54,float(page.mediabox.height)/72*2.54],'fonts':fonts,'svg_images':0,'svg_text_count':len(tree.xpath('//s:text',namespaces=ns)),'svg_tspans':math_nodes,'PNG':{'mode':png.mode,'pixels':png.size,'dpi':list(map(float,png.info['dpi']))},'TIFF':{'mode':tif.mode,'pixels':tif.size,'dpi':list(map(float,tif.info['dpi'])),'icc_bytes':len(tif.info['icc_profile'])},'pdf':'CMYK vector, embedded ICC, no images/shadings','cmyk_pdf_rendered':True}
 audit[pdf.stem]=data
assert len(audit)==3,'Expected exactly three final figure PDFs'
(r/'vector_font_and_raster_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'figures_checked':len(audit),'status':'PASS','q_chi_Arial_italic':True,'CMYK_PDF_TIFF':True,'PDF_width_cm':15},ensure_ascii=False))
