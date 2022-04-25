//
// Created by jins on 10/16/20.
//

#ifndef MAXINF_BRGSTKNN_2020_SOURCE_CONFIG_H

#define scale_coordinate 1000000.0 //用于将路网顶点坐标数值(有些坐标数值为小数)化整
#define scale_edge 1.0 //1000000.0
#define DATASET_HOME  "../../../data/"


//map， 对应数据集的选择
//#define  LV_Source1
#define  LasVegas_Source1
#define  Brightkite_Source1
#define  Gowalla_Source1
#define  Twitter_Source



#ifdef LV_Source
#define dataset_name "LV"
#define FILE_SOURCE  "../../../data/LV/source_dataset/"
//#define ROAD_PATH  "../../../data/LV/road/"
#endif


#ifdef LasVegas_Source
#define dataset_name "LasVegas"
#define FILE_SOURCE  "../../../data/LasVegas/source_dataset/"
#endif


#ifdef Brightkite_Source
#define dataset_name "Brightkite"
#define FILE_SOURCE  "../../../data/Brightkite/source_dataset/"
#endif

#ifdef Gowalla_Source
#define dataset_name "Gowalla"
#define Large_DATASET
#define FILE_SOURCE  "../../../data/Gowalla/source_dataset/"
#endif

#ifdef Twitter_Source
#define dataset_name "Twitter"
#define Large_DATASET
#define FILE_SOURCE  "../../../data/Twitter/source_dataset/"
#endif



#define MAXINF_BRGSTKNN_2020_SOURCE_CONFIG_H

#endif //MAXINF_BRGSTKNN_2020_SOURCE_CONFIG_H
