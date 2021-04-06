# SangerAlignR

#### 介绍
遗传性疾病sanger测序验证便捷查看

修改自：https://github.com/jonathonthill/sangerseqR

#### 依赖包
```r
install.packages('tinytex')
tinytex::install_tinytex()
# to uninstall TinyTeX, run tinytex::uninstall_tinytex() 
###  以下代码来自解螺旋

### 直接设定镜像 ----
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options("BioC_mirror" = "http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

#### 安装R包的一个自定义函数pkgs_in() ----
pkgs_in <- function(pkgs) {
  # 设置为国内清华镜像
  options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options("BioC_mirror" = "http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
  
  # 首先安装BiocManager，若已安装则跳过
  if (!requireNamespace("BiocManager", quietly = TRUE)) { 
    install.packages("BiocManager",ask = F, update = F)
  }
  
  # 安装stringr，若已安装则跳过
  if (!requireNamespace("stringr", quietly = TRUE)) { 
    install.packages("stringr",ask = F, update = F)
  }
  
  # 去重，识别有无github格式的安装包
  pkgs <- unique(pkgs)
  pkgs2 <- pkgs
  logi <- stringr::str_detect(pkgs2, "/")
  pkgs2[logi] <- stringr::str_match(pkgs2[logi], ".*/(.*)$")[,2]
  
  # 安装pkgs中尚未安装的包
  new <- !(sapply(pkgs2, requireNamespace, quietly = T))
  
  # 显示需安装的包
  if (sum(new) > 0) {
    cat("pkgs to install: ", pkgs[new], "\n")
  } else {
    cat("All pkgs already installed \n")
  }
  
  # install pkgs
  if(any(new)) BiocManager::install(pkgs[new], ask = F, update = F)
}

#### 3 需要安装的pkgs ----
pkgs <- c("tidyverse", "shiny", "knitr", "Biostrings", "msa", "msaR", "sangerseqR",
          "devtools", "DECIPHER", "RUnit", "tinytex")
tinytex::install_tinytex()

#### 4 安装pkgs中的R包 ----
###  pkgs_in()函数的特点 
# 可随时添加一个或多个安装包，若已安装的则不会再重复安装
# 支持安装github包，输入格式如“tidyverse/dplyr”，中间包含 / 
# 确认输入的R包名字是正确的，注意区分字母大小写
# 可反复运行，安装不成功时重启R session,再次安装
# 重启R session快捷键：Crtl + Shift + F10
# 反复运行，依然安装不成功，再次检查输入是否正确；也可能由于网络不佳，可换用其他安装方法

pkgs_in(pkgs)  
pkgs_in(pkgs)
pkgs_in(pkgs)
```

#### 直接安装修改后打包好的R包

SangerAlignR.tar.gz

```r
devtools::install_local(path="./SangerAlignR.tar.gz")
SangerAlignR::PolyPeakParser()
```

