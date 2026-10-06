"""Rebuild sex-analysis metadata from the input table, not archived outputs."""
from common import *
OUT=RUN/'sex_results';OUT.mkdir(parents=True,exist_ok=True)
m=pd.read_csv(RAW/'sample_information_10317_matched_4168.tsv',sep='\t',dtype=str).set_index('sample_name').loc[X.index]
summary={}
for c in ['contraceptive','pregnant','sex','age_years','diet_type','collection_timestamp']:
    if c in m:summary[c]={'counts':m[c].fillna('<NA>').value_counts().head(15).to_dict(),'by_sex':pd.crosstab(m['sex'],m[c].fillna('<NA>')).to_dict() if c in ['contraceptive','pregnant'] else None}
(OUT/'M01_额外协变量可用性.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
m[['contraceptive','pregnant','sex']].to_csv(OUT/'_additional_metadata.tsv',sep='\t',index_label='sample_id')
print('Sex-related metadata regenerated',len(m))
