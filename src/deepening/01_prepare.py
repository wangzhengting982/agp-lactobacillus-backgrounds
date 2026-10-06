from common import *
F=pd.read_pickle(D/'_food_primary.pkl');N=pd.read_pickle(D/'_nutrient_primary.pkl');E=N.Energy_in_kcal
groups={
'酸奶':['yogurt_all_types_except_frozen'],
'奶酪':['all_other_cheese_such_as_american_cheddar_or_cream_cheese_including_che_19c3afa9','cottage_cheese_and_ricotta_cheese','low_or_reduced_fat_cheese_including_cheese_used_in_cooking'],
'饮用牛奶':[f'milk_{v}_milk_as_a_beverage' for v in ['1','2','skim','whole']],
'水果':['all_other_fruits','apples_applesauce_and_pears','apricots_dried','apricots_fresh_or_canned','bananas','berries_such_as_strawberries_and_blueberries','cantaloupe_melon_and_mango_in_season','cherries_fresh','dried_fruit_other_than_apricots_such_as_raisins_and_prunes','grapes_fresh','oranges_grapefruit_and_tangerines_not_juice','peaches_nectarines_and_plums','pineapple_fresh_and_canned','watermelon_and_red_melon'],
'蔬菜':['broccoli','carrots_cooked','carrots_raw','cauliflower_cabbage_and_brussels_sprouts','cooked_greens_such_as_kale_mustard_greens_and_collards','cooked_greens_such_as_spinach_swiss_chard_and_beet_greens','fresh_tomatoes','green_or_string_beans','green_peas','green_peppers_and_green_chilies_cooked','green_peppers_and_green_chilies_raw','green_salad_lettuce_or_spinach','onions_and_leeks','red_peppers_and_red_chilies_cooked','red_peppers_and_red_chilies_raw','summer_squash_and_zucchini','winter_squash_such_as_acorn_butternut_and_pumpkin'],
'豆类':['all_other_beans_such_as_baked_beans_lima_beans_and_chili_without_meat','bean_soups_such_as_pea_lentil_and_black_bean','cooked_soybeans_or_edamame','refried_beans','tofu','tempeh'],
'全谷食品':['complete_or_primarily_whole_grain_cold_cereal','cooked_whole_grain_cereals','lowfat_whole_grain_crackers','regular_whole_grain_crackers','whole_grain_breads_including_bagels_and_rolls','whole_grain_breads_including_bagels_and_rolls_100_whole_grains','whole_kernel_grains_such_as_brown_rice']}
features=pd.DataFrame(index=X.index);ledger=[]
for g,names in groups.items():
 cols=[n+'__dose_g_day' for n in names];assert set(cols)<=set(F.columns)
 a=F[cols].copy();dups=a.T.duplicated();keep=a.columns[~dups]
 features[g]=scale(np.log1p(a[keep].sum(axis=1)/E*1000))
 for c in cols:ledger.append(dict(group=g,source=c,included=bool(c in keep),unit='g/day',processing='sum included item amounts; g/1000 kcal; log1p and standardize'))
for g,c in [('植物蛋白','Vegetable_Protein_in_g'),('动物蛋白','Animal_Protein_in_g'),('添加糖','Added_Sugars__by_Total_Sugars__in_g')]:
 features[g]=scale(np.log1p(N[c]/E*1000));ledger.append(dict(group=g,source=c,included=True,unit='g/day',processing='g/1000 kcal; log1p and standardize'))
for g,c in [('植物发酵食品','plant_ferment'),('益生菌使用','probiotic')]:features[g]=X[c];ledger.append(dict(group=g,source=c,included=True,unit='ordinal 0 to 4',processing='same standardized original ordinal'))
assert features.shape==(2748,12) and np.isfinite(features).all().all()
features.to_pickle(O/'_diet12.pkl');save(ledger,'D01_12项摄入定义');features.describe().T.to_csv(O/'D01_摄入变量描述.tsv',sep='\t')
tax=pd.read_csv(RAW/'4168_annotated_feature_table.tsv',sep='\t',index_col=0)
forbidden=tax.index.str.contains(r'(?:^|;)o__Lactobacillales(?:;|$)|(?:^|;)f__Bifidobacteriaceae(?:;|$)');ok=tax.index.str.startswith('d__Bacteria;')&~forbidden
ge=tax.index.to_series().str.extract(r'(?:^|;)g__([^;]+)',expand=False);fa=tax.index.to_series().str.extract(r'(?:^|;)f__([^;]+)',expand=False).fillna('unclassified_family');known=ge.notna()&~ge.fillna('').str.contains('unclassified|uncultured',case=False)
allbg=tax.loc[ok&known].groupby((fa+'|'+ge)[ok&known]).sum().T;allbg.columns=['microbe::'+c for c in allbg]
den=allbg.sum(axis=1);lg=np.log(allbg.reindex(columns=Z.columns).div(den,axis=0)+.00005);clr=lg.sub(lg.mean(axis=1),axis=0);allz=clr.sub(clr.loc[X.index].mean()).div(clr.loc[X.index].std(ddof=0))
assert np.allclose(allz.loc[X.index,Z.columns],Z,atol=1e-10)
fm=[c for c in allbg if c.startswith('microbe::Peptoniphilaceae|')];frac=allbg[fm].sum(axis=1)/den;fv=np.log((frac+.00005)/(1-frac+.00005));fv=(fv-fv.loc[X.index].mean())/fv.loc[X.index].std(ddof=0)
pd.DataFrame({'family_balance':fv,'named_depth':den,'depth_nonlab':np.log10(tax.loc[ok].sum(axis=0)+1)}).to_pickle(O/'_all_depth_family.pkl');allz[BC].to_pickle(O/'_all_core_clr.pkl')
(O/'准备核验.json').write_text(json.dumps({'diet_features':features.shape[1],'CLR_max_error':float(np.abs(allz.loc[X.index,Z.columns]-Z).to_numpy().max()),'family_members':fm,'family_positive_n':int(allbg.loc[X.index,fm].sum(axis=1).gt(0).sum())},ensure_ascii=False,indent=2),encoding='utf-8')
print('prepared',features.shape,'family',len(fm),flush=True)
