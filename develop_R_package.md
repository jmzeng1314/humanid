# 开发R包的经验分享
首先需要安装一些必备的包，然后需要抽空学习 Hadley Wickham 的 [R package ](http://r-pkgs.had.co.nz/r.html) 这本书
```r
install.packages('devtools')
## 检查R代码是否规范 
install.packages("formatR")
formatR::tidy_dir("R")
## 检查R包写的怎么样
install.packages("lintr")
lintr::lint_package()
```

## 用BiocManager来管理和安装github上面的包

BiocManager 负责解决biocondutor的依赖 然后还能调用devtools 和install.packages,这样写在R包里面的依赖包关系会被自动下载及更新。
```r
chooseCRANmirror()
install.packages("BiocManager")
BiocManager::install("jmzeng1314/humanid")
```

这里我以 https://github.com/jmzeng1314/humanid 为例子来讲解

### 如何与GitHub同步

首先是可以使用Rstudio编辑器的同步功能，需要自行摸索设置。

然后是可以使用Git命令行来操作，需要背诵的命令行有点多，比如Rstudio编辑器的同步功能方便。

### R包的文件夹结构

因为Rstudio编辑器自带新建R包功能，这里就不演示从头开始创建R包了，Rstudio编辑器的新建R包就会自动化把R包的文件格式都设置好。最重要的就3个文件夹，分别是R,man,data，接下来我们就一一介绍：

#### 首先是R代码文件夹

#### 然后是man文件夹

#### 接着是data文件夹

### 关于DESCRIPTION

### 关于NAMESPACE

