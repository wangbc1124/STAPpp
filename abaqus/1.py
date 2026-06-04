# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *

# 1. 打开你的 ODB 文件（请将下面的路径和文件名替换为你的实际文件）
odb_path = 'Bridge-1.odb'  # 替换为你的 ODB 文件路径
odb = session.openOdb(name=odb_path)

# 获取装配体下的所有实体列表
all_instances = odb.rootAssembly.instances

# 2. 创建一个总的 CSV 文件
output_file = 'all_instances_nodes.csv'
with open(output_file, 'w') as f:
    # 写入表头，增加了一列 实体名称(Instance_Name)
    f.write('Instance_Name, Node_ID, X, Y, Z\n')
    
    # 3. 循环遍历每一个实体
    for inst_name in all_instances.keys():
        instance = all_instances[inst_name]
        nodes = instance.nodes
        
        # 逐个写入该实体的节点坐标
        for node in nodes:
            node_id = node.label
            x, y, z = node.coordinates # 初始坐标（如需变形后坐标需加位移U）
            f.write('%s, %d, %f, %f, %f\n' % (inst_name, node_id, x, y, z))

print('所有实体的节点坐标已成功分类导出至: %s' % output_file)