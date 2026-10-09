当前三张主图复现（图3合并检出状态与检出者丰度）

运行：python rebuild_figures.py --output-dir ../combined_figure3
依赖：requirements.txt；Windows 默认 Arial、幼圆字体。
可通过 --source-dir、--font-latin、--font-cjk 指定输入与字体路径。

输入 source_data 原字节复制自上一轮已验证当前图件包；未拟合模型或修改数据。
图1、图2沿用上一轮排版；图3为15×11 cm，上方(a)为2748人四种检出状态，
下方(b)为502名乳酪杆菌属检出者及469名乳杆菌属检出者的丰度关联。
默认只生成三张主图的可编辑SVG、600 dpi RGB PNG、600 dpi CMYK TIFF及预览，
不生成S1或PDF。图1/2沿用全部数值和轴，图3沿用原S1与图3的全部数值及轴。
所有66行绘图值须与original_plotted_values.json精确一致；39条CI与实际绘图对象核对。
输出numeric_plot_audit.json、typography_audit.json、export_audit.json记录数值、字体和导出核验。
所有文字7 pt、黑色，Arial与幼圆；PNG/TIFF色彩转换使用随包ICC。

当前30张结果表的对应关系：冻结绘图JSON保留整合前的工作表名与行号；其中C13的两目标效应量现并入C12，A01原始检出者结果位于合并表的all_detected方案，S02概率列不变。已按模型键逐项核对当前工作簿，54行统计绘图值完全一致，另12行图1来源数据保持不变。旧行号用于追踪来源，不是当前Excel的单元格地址。
