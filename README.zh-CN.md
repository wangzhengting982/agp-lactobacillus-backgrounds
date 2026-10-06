# AGP乳酸菌属肠道菌群背景分析

**《乳杆菌属与乳酪杆菌属的肠道菌群背景差异》配套分析代码。**

[English](README.md) · [数据与输入要求](DATA_AVAILABILITY.md) · [核验结果](reports/final_validation.json) · [结果表与代码对应关系](reports/table_map.csv)

## 项目简介

本研究使用美国肠道计划（American Gut Project）的饮食和粪便微生物组资料，比较6个乳酸菌属所关联的肠道菌群背景，重点分析乳杆菌属与乳酪杆菌属之间的差异。代码包括人群筛选、背景菌关联模型、饮食调整、敏感性分析及主图生成。

| 项目 | 内容 |
|---|---|
| 主分析人群 | 2748名成年人 |
| 背景菌属 | 156个 |
| 目标菌属 | 6个 |
| 已核验输出 | 40张结果表、3幅主图 |
| 已验证环境 | Windows、Python 3.12 |

## 代码结构

- [`src/preparation/`](src/preparation/)：筛选人群、构建模型输入。
- [`src/core/`](src/core/)：属间差异及分类敏感性分析。
- [`src/deepening/`](src/deepening/)：饮食主成分、饮食调整、检出状态、丰度等分析。
- [`src/sensitivity/`](src/sensitivity/)：其他稳健性分析。
- [`src/figures/`](src/figures/)：绘图程序和汇总绘图参考值。
- [`scripts/`](scripts/)：完整运行入口及核验工具。
- [`reports/`](reports/)：本地运行记录、核验结果和输入清单。

## 如何运行

先按照[输入说明](DATA_AVAILABILITY.md)准备 `data/archive/` 下的完整输入包和 `reference/` 下的3个参考文件，再在项目目录中运行：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
.\.venv\Scripts\python.exe scripts/run_all.py
```

需要Python 3.12。程序依次检查文件哈希、重建人群与模型输入、拟合模型、核对40张结果表并重绘3幅主图。分析输出保存到 `runs/`，运行与核验记录保存到 `reports/`。已完成环境配置的本地项目也可双击 `run_local.cmd`。

## 已完成的验证

2026年10月6日的本地运行通过了40张当前稿件结果表的核对，最大绝对数值差为8.24×10⁻¹³。3幅主图均由新计算的模型结果重新绘制，按200 dpi渲染后与投稿PDF逐像素一致。

详细证据见[核验汇总](reports/final_validation.json)、[完整运行记录](reports/full_run_status.json)和[中文核验说明](docs/本地复现验证报告.md)。展示用报告中的本机解释器路径已替换为通用占位符，计算结果、版本、时间及退出码保持原记录。

## 数据与复现范围

仓库提供代码、汇总绘图参考值及核验报告。个体级数据、完整输入包、正文、补充材料和参考工作簿另行保存；运行完整流程需同时准备这些输入文件。

本流程从处理后的表格、BIOM计数、已有分类注释和固定纵向样本索引开始，重新计算人群、模型、结果表和图片。原始测序数据处理、最初问卷匹配及序列分类位于本流程上游。研究结论为横断面关联。
