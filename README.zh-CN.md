# PolypAtlas — CRSwNP 多模态转录组资源

[English](README.md) | **简体中文**

**在线资源：https://polyp.erbihou.org/**

PolypAtlas 是一个为慢性鼻窦炎伴鼻息肉（CRSwNP）策展、统一协调的多模态转录组资源，完全基于公开数据构建。它整合了 GEO 上目录的全部人类 CRS / 鼻息肉公开转录组数据集（策展 66 个系列，核心数据集 32 个），涵盖 bulk 转录组、单细胞 RNA 测序与空间转录组，并在同一套统一流程下重新处理，使任何基因都能跨队列、跨细胞类型、跨组织区室进行比较。

## 内容

- **Bulk 转录组：** 16 个队列、30 个对比、19,266 个基因
- **单细胞 RNA 测序：** 10 个数据集、476,254 个细胞、105 位供体
- **空间转录组（GeoMx）：** 息肉与对照的组织区室
- **内置验证：** marker 方向回收率 87.5%；8 个 marker 中 7 个被 dupilumab 逆转

## 网站功能

- **基因查询** — 任一基因在 16 个队列中的 bulk logFC/FDR、单细胞对比与空间区室数值，附通路成员与共识注释
- **通路查询** — 检索一条富集通路，查看它出现在哪些队列与细胞类型中
- **队列浏览器** — 每个数据集的完整元数据，每个数据集的局限就写在其数据旁边
- **下载** — 证据矩阵、差异表达表、组成表与策展表
- **方法** — 完整流程与已知局限，全文写明

## 数据可用性

全部处理后数据存档于 Zenodo（CC BY 4.0）：**DOI [10.5281/zenodo.23448779](https://doi.org/10.5281/zenodo.23448779)**。同样的表格也可直接从网站下载页获取。

## 引用

如在工作中使用 PolypAtlas，请引用：

> PolypAtlas: a curated multi-modal transcriptomic resource for chronic rhinosinusitis with nasal polyps. https://polyp.erbihou.org/ — archived data: Zenodo, DOI 10.5281/zenodo.23448779.

## 仓库结构

- `build/` — 静态站点生成器
- `site/` — 生成的网站（经 `gh-pages` 分支在上述地址提供服务）

## 联系方式

Peiliang Lei（of@live.com），石家庄市第三医院耳鼻喉科，中国河北石家庄。
