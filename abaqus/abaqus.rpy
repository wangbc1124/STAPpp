# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2024 replay file
# Internal Version: 2023_09_21-20.55.25 RELr426 190762
# Run by Nie_J on Wed Jun 10 15:22:26 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=306.0791015625, 
    height=125.124992370605)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='Bridge-1.odb')
#: 模型: D:/STAPpp-master/STAPpp-master/abaqus/Bridge-1.odb
#: 装配件个数:         1
#: 装配件实例个数: 0
#: 部件实例的个数:     27
#: 网格数:             27
#: 单元集合数:       10
#: 结点集合数:          104
#: 分析步的个数:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1256.65, 
    farPlane=1541.35, width=779.121, height=344.785, cameraPosition=(-86.2781, 
    -1278.92, 606.546), cameraUpVector=(-0.171556, 0.694467, 0.698773), 
    cameraTarget=(23.9656, -25.1282, 44.9946))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1086.94, 
    farPlane=1680.05, width=673.902, height=298.222, cameraPosition=(-666.277, 
    -744.544, 1002.63), cameraUpVector=(0.211843, 0.90755, 0.362595), 
    cameraTarget=(15.3491, -17.1895, 50.8789))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1144.86, 
    farPlane=1622.04, width=709.815, height=314.115, cameraPosition=(-339.379, 
    -287.541, 1355.58), cameraUpVector=(0.252383, 0.964123, -0.0822758), 
    cameraTarget=(16.5955, -15.447, 52.2247))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1153.47, 
    farPlane=1613.44, width=609.279, height=269.624, viewOffsetX=-7.83755, 
    viewOffsetY=-2.82443)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1183.2, 
    farPlane=1587.99, width=624.984, height=276.574, cameraPosition=(-190.413, 
    -397.76, 1359.12), cameraUpVector=(0.151076, 0.987781, -0.0382644), 
    cameraTarget=(17.2332, -15.0818, 51.4824), viewOffsetX=-8.03958, 
    viewOffsetY=-2.89723)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1302.18, 
    farPlane=1496.93, width=687.829, height=304.385, cameraPosition=(-30.9959, 
    -1383.68, 254.004), cameraUpVector=(-0.197359, 0.473506, 0.858395), 
    cameraTarget=(16.8957, -21.6579, 48.8925), viewOffsetX=-8.848, 
    viewOffsetY=-3.18856)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
o1 = session.openOdb(name='D:/STAPpp-master/STAPpp-master/abaqus/Bridge-2.odb', 
    readOnly=False)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: 模型: D:/STAPpp-master/STAPpp-master/abaqus/Bridge-2.odb
#: 装配件个数:         1
#: 装配件实例个数: 0
#: 部件实例的个数:     27
#: 网格数:             27
#: 单元集合数:       10
#: 结点集合数:          104
#: 分析步的个数:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1229.42, 
    farPlane=1591.99, width=762.234, height=337.312, cameraPosition=(481.692, 
    -1323.18, -41.6112), cameraUpVector=(-0.147831, 0.23902, 0.959696), 
    cameraTarget=(24.9958, -25.6215, 43.5719))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1217.12, 
    farPlane=1604.83, width=754.61, height=333.938, cameraPosition=(494.579, 
    -1313.71, 188.628), cameraUpVector=(-0.252749, 0.361869, 0.897312), 
    cameraTarget=(25.2926, -25.4033, 48.8742))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1279.81, 
    farPlane=1533.34, width=793.475, height=351.137, cameraPosition=(91.4213, 
    -1349.35, -341.791), cameraUpVector=(-0.046088, 0.0594423, 0.997167), 
    cameraTarget=(15.9324, -26.2308, 36.5593))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1305.12, 
    farPlane=1505.35, width=809.167, height=358.081, cameraPosition=(-15.5216, 
    -1400.35, 161.806), cameraUpVector=(0.0830987, 0.408979, 0.908752), 
    cameraTarget=(13.777, -27.2588, 46.7094))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1212.83, 
    farPlane=1608.92, width=751.948, height=332.76, cameraPosition=(534.453, 
    -1305.48, 16.2825), cameraUpVector=(-0.0640794, 0.312908, 0.947619), 
    cameraTarget=(24.3497, -25.435, 43.9119))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1041.64, 
    farPlane=1799.1, width=645.812, height=285.792, cameraPosition=(1341.5, 
    -89.905, 508.113), cameraUpVector=(-0.608791, 0.19041, 0.770141), 
    cameraTarget=(43.0863, 2.7862, 55.3304))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1165.75, 
    farPlane=1668.4, width=722.759, height=319.843, cameraPosition=(710.471, 
    -1149.83, 470.653), cameraUpVector=(-0.0243327, 0.670893, 0.741155), 
    cameraTarget=(23.6391, -29.8788, 54.176))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1175.12, 
    farPlane=1659.65, width=728.567, height=322.413, cameraPosition=(712.396, 
    -1206.39, 259.667), cameraUpVector=(0.632701, 0.77224, 0.057752), 
    cameraTarget=(23.6919, -31.4297, 48.3905))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1282.36, 
    farPlane=1547.46, width=795.057, height=351.837, cameraPosition=(86.1321, 
    -1325.81, 531.544), cameraUpVector=(0.907797, 0.306385, -0.286415), 
    cameraTarget=(6.3867, -34.7296, 55.9031))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1105.15, 
    farPlane=1717.61, width=685.188, height=303.216, cameraPosition=(-691.71, 
    -96.9908, 1271.35), cameraUpVector=(0.69073, -0.721426, -0.0493626), 
    cameraTarget=(-13.786, -2.86116, 75.0894))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1235.9, 
    farPlane=1586.6, width=766.25, height=339.089, cameraPosition=(-209.279, 
    1171.74, 803.114), cameraUpVector=(-0.0394175, -0.797267, 0.602338), 
    cameraTarget=(-2.44984, 26.9515, 64.0868))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1273.91, 
    farPlane=1548.59, width=213.565, height=94.5089, viewOffsetX=45.5194, 
    viewOffsetY=-32.9793)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1103.67, 
    farPlane=1776.28, width=185.024, height=81.8789, cameraPosition=(-1370.89, 
    -57.2654, 486.051), cameraUpVector=(0.615414, -0.155034, 0.772807), 
    cameraTarget=(-64.2039, 55.5099, 62.5767), viewOffsetX=39.4363, 
    viewOffsetY=-28.572)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1134.87, 
    farPlane=1742.95, width=190.255, height=84.1937, cameraPosition=(-1089.74, 
    658.197, 716.19), cameraUpVector=(0.310799, -0.780904, 0.541843), 
    cameraTarget=(-29.2794, 50.2326, 79.5728), viewOffsetX=40.5512, 
    viewOffsetY=-29.3798)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1105.7, 
    farPlane=1772.12, width=574.4, height=254.19, viewOffsetX=-19.7982, 
    viewOffsetY=-28.0142)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1271.6, 
    farPlane=1680.72, width=660.583, height=292.328, cameraPosition=(-409.073, 
    1184, 826.073), cameraUpVector=(-0.253068, -0.837587, 0.484154), 
    cameraTarget=(-21.5995, 89.3475, 83.7323), viewOffsetX=-22.7687, 
    viewOffsetY=-32.2174)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1351.27, 
    farPlane=1647.86, width=701.968, height=310.642, cameraPosition=(154.005, 
    1359.12, 660.029), cameraUpVector=(-0.566805, -0.580931, 0.584167), 
    cameraTarget=(10.5618, 119.657, 74.6781), viewOffsetX=-24.1952, 
    viewOffsetY=-34.2358)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1242.07, 
    farPlane=1724.15, width=645.241, height=285.539, cameraPosition=(597.525, 
    -1067.48, 886.199), cameraUpVector=(-0.388322, 0.708791, 0.588915), 
    cameraTarget=(98.3685, -45.7644, 107.464), viewOffsetX=-22.2399, 
    viewOffsetY=-31.4691)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1309.15, 
    farPlane=1621.72, width=680.091, height=300.961, cameraPosition=(160.194, 
    -1065.91, 1040.97), cameraUpVector=(0.0828368, 0.882743, 0.462496), 
    cameraTarget=(83.689, -51.2639, 111.384), viewOffsetX=-23.4411, 
    viewOffsetY=-33.1687)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1286.98, 
    farPlane=1633.56, width=668.573, height=295.864, cameraPosition=(-342.169, 
    -1406.74, 251.039), cameraUpVector=(0.0945251, 0.440805, 0.892612), 
    cameraTarget=(49.1978, -96.7182, 77.3335), viewOffsetX=-23.0441, 
    viewOffsetY=-32.6069)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1321.65, 
    farPlane=1608.94, width=686.583, height=303.834, cameraPosition=(-205.706, 
    -1439.85, 238.356), cameraUpVector=(0.0111361, 0.446777, 0.894576), 
    cameraTarget=(55.4733, -96.1836, 77.585), viewOffsetX=-23.6649, 
    viewOffsetY=-33.4852)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1152.93, 
    farPlane=1872.47, width=598.933, height=265.046, cameraPosition=(1505.4, 
    149.956, 27.7222), cameraUpVector=(-0.281062, -0.267008, 0.921798), 
    cameraTarget=(134.57, 13.3602, 68.5614), viewOffsetX=-20.6438, 
    viewOffsetY=-29.2104)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-1', frame=1)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1164.24, 
    farPlane=1854.42, width=604.809, height=267.647, cameraPosition=(1428.52, 
    485.348, -6.08458), cameraUpVector=(-0.197378, -0.306893, 0.931053), 
    cameraTarget=(125.018, 43.4334, 65.3434), viewOffsetX=-20.8463, 
    viewOffsetY=-29.497)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-1', frame=0)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-1', frame=1)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1215.78, 
    farPlane=1760.9, width=631.583, height=279.495, cameraPosition=(-757.731, 
    1040.13, 793.622), cameraUpVector=(0.974219, 0.124517, 0.188132), 
    cameraTarget=(-34.1326, 74.6022, 127.535), viewOffsetX=-21.7691, 
    viewOffsetY=-30.8027)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1123.43, 
    farPlane=1862.24, width=583.607, height=258.264, cameraPosition=(-1427.37, 
    8.53048, 483.785), cameraUpVector=(0.536254, 0.516697, 0.667425), 
    cameraTarget=(-99.9203, -1.44307, 113.269), viewOffsetX=-20.1155, 
    viewOffsetY=-28.4629)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1183.64, 
    farPlane=1797.49, width=614.884, height=272.105, cameraPosition=(-927.784, 
    -785.812, 908.095), cameraUpVector=(0.0354724, 0.935704, 0.350998), 
    cameraTarget=(-60.8641, -39.3785, 139.465), viewOffsetX=-21.1935, 
    viewOffsetY=-29.9883)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-1', frame=0)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1168.35, 
    farPlane=1808.66, width=606.943, height=268.591, cameraPosition=(-1032.9, 
    -742.016, 824.238), cameraUpVector=(0.0543624, 0.930881, 0.361254), 
    cameraTarget=(-70.3596, -36.982, 134.353), viewOffsetX=-20.9198, 
    viewOffsetY=-29.601)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-1', frame=0)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    deformationScaling=UNIFORM, uniformScaleFactor=10)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1129.51, 
    farPlane=1853.56, width=586.768, height=259.663, cameraPosition=(-1427.02, 
    -272.662, 389.48), cameraUpVector=(0.261611, 0.869828, 0.418281), 
    cameraTarget=(-105.307, -5.14802, 104.841), viewOffsetX=-20.2244, 
    viewOffsetY=-28.617)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-1', frame=1)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1326.96, 
    farPlane=1647.59, width=689.34, height=305.054, cameraPosition=(-173.476, 
    -1130.21, 1001.83), cameraUpVector=(-0.506397, 0.78703, 0.352343), 
    cameraTarget=(-14.9948, -61.4514, 146.172), viewOffsetX=-23.7598, 
    viewOffsetY=-33.6195)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1232.02, 
    farPlane=1759.56, width=640.022, height=283.229, cameraPosition=(-751.9, 
    -1167.97, -506.451), cameraUpVector=(-0.157781, 0.0316171, 0.986968), 
    cameraTarget=(-49.8738, -114.022, 37.4868), viewOffsetX=-22.0599, 
    viewOffsetY=-31.2142)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1317.74, 
    farPlane=1677.08, width=684.552, height=302.935, cameraPosition=(-238.489, 
    -1194.79, -821.547), cameraUpVector=(-0.207591, -0.255863, 0.944161), 
    cameraTarget=(-3.3043, -119.667, 8.09334), viewOffsetX=-23.5947, 
    viewOffsetY=-33.386)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1325.32, 
    farPlane=1665.69, width=688.489, height=304.678, cameraPosition=(-339.464, 
    -1428.05, -239.712), cameraUpVector=(0.099575, 0.102074, 0.989781), 
    cameraTarget=(-0.291209, -125.342, 55.9724), viewOffsetX=-23.7304, 
    viewOffsetY=-33.578)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1399.63, 
    farPlane=1590.37, width=727.094, height=321.761, cameraPosition=(93.6555, 
    -1488.07, 167.871), cameraUpVector=(-0.0670441, 0.382129, 0.921674), 
    cameraTarget=(31.7961, -113.325, 91.8498), viewOffsetX=-25.061, 
    viewOffsetY=-35.4607)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1295.01, 
    farPlane=1694.41, width=672.745, height=297.71, cameraPosition=(-537.109, 
    -1392.36, -44.7451), cameraUpVector=(0.138969, 0.215317, 0.966605), 
    cameraTarget=(-16.8214, -121.518, 72.6467), viewOffsetX=-23.1877, 
    viewOffsetY=-32.8101)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1120.38, 
    farPlane=1863.35, width=582.029, height=257.566, cameraPosition=(-1425.72, 
    -48.4044, 488.995), cameraUpVector=(0.569761, -0.317737, 0.757902), 
    cameraTarget=(-102.345, -39.202, 104.149), viewOffsetX=-20.061, 
    viewOffsetY=-28.3858)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1204.24, 
    farPlane=1777.05, width=625.597, height=276.846, cameraPosition=(-909.459, 
    938.379, 768.716), cameraUpVector=(0.186046, -0.77545, 0.603377), 
    cameraTarget=(-90.3703, 38.6899, 121.286), viewOffsetX=-21.5626, 
    viewOffsetY=-30.5106)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1335.3, 
    farPlane=1641.2, width=693.683, height=306.976, cameraPosition=(-24.3349, 
    46.7924, 1538.11), cameraUpVector=(0.176915, -0.944351, -0.277313), 
    cameraTarget=(-21.0464, -35.0639, 162.318), viewOffsetX=-23.9093, 
    viewOffsetY=-33.8311)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1380.86, 
    farPlane=1595.64, width=64.2271, height=28.4224, viewOffsetX=-23.2169, 
    viewOffsetY=-28.0554)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1378.68, 
    farPlane=1597.26, width=64.1258, height=28.3776, cameraPosition=(-32.8353, 
    105.881, 1534.63), cameraUpVector=(0.184536, -0.953449, -0.238499), 
    cameraTarget=(-21.4787, -30.7131, 163.236), viewOffsetX=-23.1803, 
    viewOffsetY=-28.0112)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1316.48, 
    farPlane=1659.46, width=775.809, height=343.319, viewOffsetX=-58.4461, 
    viewOffsetY=19.6765)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1319.06, 
    farPlane=1570.04, width=777.326, height=343.991, cameraPosition=(186.32, 
    1432.62, 46.7419), cameraUpVector=(-0.0514261, -0.338215, 0.939663), 
    cameraTarget=(-1.79015, 67.3478, 33.9487), viewOffsetX=-58.5605, 
    viewOffsetY=19.715)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1134.23, 
    farPlane=1834.25, width=668.404, height=295.789, cameraPosition=(1124.05, 
    640.465, 778.863), cameraUpVector=(0.115986, -0.992971, 0.023566), 
    cameraTarget=(62.9465, 38.6392, 137.467), viewOffsetX=-50.3547, 
    viewOffsetY=16.9524)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1230.71, 
    farPlane=1725.83, width=725.257, height=320.948, cameraPosition=(444.785, 
    -79.2974, 1459.12), cameraUpVector=(0.498222, -0.681762, -0.535701), 
    cameraTarget=(-2.45594, -55.7997, 155.679), viewOffsetX=-54.6378, 
    viewOffsetY=18.3943)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1320.3, 
    farPlane=1584.99, width=778.051, height=344.311, cameraPosition=(31.3809, 
    1301.89, 695.779), cameraUpVector=(0.728748, -0.574448, 0.372741), 
    cameraTarget=(-28.5366, 51.7032, 118.727), viewOffsetX=-58.6151, 
    viewOffsetY=19.7333)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1282.79, 
    farPlane=1622.51, width=1250.3, height=553.295, viewOffsetX=-118.976, 
    viewOffsetY=52.2203)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1277.5, 
    farPlane=1627.79, width=1245.15, height=551.016, cameraPosition=(13.6959, 
    1332.56, 631.159), cameraUpVector=(0.35501, -0.680406, 0.641105), 
    cameraTarget=(-46.2216, 82.3777, 54.1067), viewOffsetX=-118.486, 
    viewOffsetY=52.0052)
