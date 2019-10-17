<!--
 * FileName: README.md
 * Version: 1.0
 * License: MIT License
* 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-17 13:55:35
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-17 14:02:20
 * Description: README
 -->

# Distribution

## 说明

此程序用于从 HepMC2 文件中挑选 `Z` Boson 以及符合条件的 `Jet`，并计算 `Z` 玻色子与 `Jet` 之间的 delta phi 角的分布，最后保存在由 `--name=[...]` 命名的指定文件内。

## 编译

1. 新建 ```build``` 目录
    
    ```mkdir build && cd build```

2. 创建 ```Makefile``` 文件
    
    ```cmake ..```

    若 `HepMC2` 和 `FastJet` 非标准安装，请使用 `-DHepMC_INCLUDE_DIRS=...` 和 `-DHepMC_LIBRARY_DIRS=...` 指定 `HepMC2` 的头文件目录和库文件目录；用 `-DFastjet_CONFIG=...` 来指定 `fastjet-config` 文件的位置。

3. 编译
    
    ```make```

## 使用

编译完成之后会生成 ```dist``` 可执行文件，然后执行
    
```./dist -n delta_phi.dat *.hepmc```

这里 `-n` 后接指定的输出文件名，`*.hepmc` 为待处理 `HepMC2` 文件。
