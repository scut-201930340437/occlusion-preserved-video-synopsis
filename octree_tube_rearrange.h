#pragma once
#include <opencv2/opencv.hpp>
#include "segment.h"
#include "defines.h"
struct OctreeNode
{
	//T data; //�ڵ�����
	int tube_id;
	std::vector<cv::Vec6d>tubes_cub;
	std::vector<OctreeNode*>tubes_oct;
	double xmin, xmax; //�ڵ����꣬��������8�����������  
	double ymin, ymax;
	double zmin, zmax;
	//OctreeNode* top_left_front, * top_left_back; //�ýڵ��8���ӽ��  
	//OctreeNode* top_right_front, * top_right_back;
	//OctreeNode* bottom_left_front, * bottom_left_back;
	//OctreeNode* bottom_right_front, * bottom_right_back;
	std::vector<OctreeNode*>children;
	OctreeNode();
	OctreeNode(double xminValue, double xmaxValue, double yminValue, double ymaxValue, double zminValue, double zmaxValue,int _tube_id=0);
};

void createOctree(OctreeNode* root,  double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
int cal(int num);
void find(OctreeNode* p, double x, double y, double z);

void insert_tube(OctreeNode* root,  OctreeNode* tube_oct);
void octree_tube_rerrange(Context& context, std::map<int, int>& tubes_group_id, double speed);