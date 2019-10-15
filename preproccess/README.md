<!--
 * FileName: README
 * Version: 1.0
 * License: MIT License
* 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-15 18:02:58
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-15 19:19:23
 * Description: README
 -->
# Preproccess

## 说明

此程序用于 HepMC2 文件的预处理，会将符合条件的事件筛选出来，并保存到程序所在目录下的 `Output` 
文件夹之内。

## 编译

1. 新建 ```build``` 目录
    
    ```mkdir build && cd build```

2. 创建 ```Makefile``` 文件
    
    ```cmake ..```

3. 编译
    
    ```make```

## 使用

编译完成之后会生成 ```prep``` 可执行文件，首先在可执行文件所在的目录新建一个 `Output` 
文件夹，否则会运行出错。然后执行
    
```mpiexec -np N ./prep *.hepmc```

这里 `N` 为 `MPI` 运行参数，`*.hepmc` 为待处理 `HepMC` 文件。
