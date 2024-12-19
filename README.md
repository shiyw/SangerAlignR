# SangerAlignR

## 介绍
遗传性疾病sanger测序验证便捷查看

Fork from https://gitee.com/x2yline/sanger-align-r

## 使用
### 依赖
可以使用renv来管理依赖，也可以手动安装

#### renv
```r
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::restore()
```

#### 手动安装
```r
install.packages("pak")
pak::pak(c("tinytex", "tidyverse", "shiny", "knitr", "Biostrings", "msa", "msaR", "sangerseqR",
          "devtools", "DECIPHER", "RUnit", "tinytex"))
tinytex::install_tinytex()
```

### 安装SangerAlignR
```r
pak::pak("shiyw/SangerAlignR")
```

### 打开ShinyAPP
```r
SangerAlignR::PolyPeakParser()
```

