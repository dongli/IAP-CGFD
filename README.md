简介
==

IAP-CGFD是中国科学院大气物理研究所王斌课题组在中国科学院大学开设的*计算地球流体力学*课程辅助文件仓库。目前包含有若干简单的计算算例供大家动手实践。

安装要求
====

演示计算程序使用了董理编写的[GEOMTK](https://github.com/dongli/geomtk)模式工具库，而GEOMTK又依赖一些C++软件库，因此需要先将依赖库安装好。另外，由于采用了`c++11`语言特性，因此需要较高版本的GCC。可以使用[PACKMAN](https://github.com/dongli/packman)来安装，比如：
```
$ packman install gcc
$ packman install armadillo boost mlpack netcdf gsl
```

成员
==

- 王斌
- 刘娟娟
- 公敬
- 李艺苑
- 董理 <dongli@lasg.iap.ac.cn>