# project1
此项目是关于学习PY，以及使用PY实现卫星定位算法，目前v1.0版本实现将rinex格式数据进行读取，并实现基于测距的单点定位算法。

## Py_study目录：
主要是关于py的一些学习，功能实现的方式等。  

## positon目录：
学习并实现卫星的定位算法

### 1、common目录：
- CoordinateTransformation.py：XYZ坐标系与BLH坐标系互相转换的功能类。  
- NAVDATA.py：包含将rinex格式的星历文件进行提取转换的方法，并打包成星历类。  
- OBSDATA.py：包含将rinex格式的观测量文件进行提取转换的方法，并打包成观测值类。  
- SATPOS.py：用于计算卫星坐标的方法。
- tool.py：包含了一些常用的常参、数学方法等工具。 
### 2、singlepoint目录
- singlepoint.py：主函数。  
- SINGLEPOINTPOSITION.py：实现单点定位算法。 
### 3、rtkpoint目录
- 还在完善

