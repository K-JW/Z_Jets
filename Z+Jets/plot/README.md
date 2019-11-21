<!--
 * FileName: README.md
 * Version: 1.0
 * License: MIT License
* 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-17 20:36:09
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-17 21:57:40
 * Description: README
 -->

# Plot

此文件夹内存放绘图代码及工具。

## 1. make-plots

`make-plots` 提取自 `Rivet-2.7.0.tar.gz`，我又增加了绘制点集的功能。这是一个 `Python` 脚本，它能读取 `*.dat` 中的数据，然后转换为 `LaTeX` 中的绘图语句，使用 `LaTeX` 引擎完成编译，输出图像。默认生成 `PDF` 文件。

[make-plots manual](https://rivet.hepforge.org/make-plots.html)