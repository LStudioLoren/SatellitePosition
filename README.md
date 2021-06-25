# project1
学习PY，以及使用PY实现卫星定位算法  
# Py_study目录：
主要是关于py的一些学习，功能实现的方式等。  

# singlepoint目录：
学习并实现卫星的基于伪距的单点定位算法实现。 

CoordinateTransformation.py：XYZ坐标系与BLH坐标系互相转换的功能类。  
NAVDATA.py：包含将rinex格式的星历文件进行提取转换的方法，并打包成星历类。  
OBSDATA.py：包含将rinex格式的观测量文件进行提取转换的方法，并打包成观测值类。  
SATPOS.py：用于计算卫星坐标的方法。  
singlepoint.py：主函数。  
SINGLEPOINTPOSITION.py：实现单点定位算法。  
tool.py：包含了一些常用的常参、数学方法等工具。 
