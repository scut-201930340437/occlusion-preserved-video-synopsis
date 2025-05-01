#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <vector>

#include "octree_tube_rearrange.h"
#include "frame_seq_fusion.h"
#include "loadata.h"
//定义八叉树节点类

//OctreeNode //节点类  
	//(T nodeValue = T(),
	//	T xminValue = T(), T xmaxValue = T(),
	//	T yminValue = T(), T ymaxValue = T(),
	//	T zminValue = T(), T zmaxValue = T(),
	//	OctreeNode<T>* top_left_front_Node = NULL,
	//	OctreeNode<T>* top_left_back_Node = NULL,
	//	OctreeNode<T>* top_right_front_Node = NULL,
	//	OctreeNode<T>* top_right_back_Node = NULL,
	//	OctreeNode<T>* bottom_left_front_Node = NULL,
	//	OctreeNode<T>* bottom_left_back_Node = NULL,
	//	OctreeNode<T>* bottom_right_front_Node = NULL,
	//	OctreeNode<T>* bottom_right_back_Node = NULL)
	//	:data(nodeValue),
	//	xmin(xminValue), xmax(xmaxValue),
	//	ymin(yminValue), ymax(ymaxValue),
	//	zmin(zminValue), zmax(zmaxValue),
	//	top_left_front(top_left_front_Node),
	//	top_left_back(top_left_back_Node),
	//	top_right_front(top_right_front_Node),
	//	top_right_back(top_right_back_Node),
	//	bottom_left_front(bottom_left_front_Node),
	//	bottom_left_back(bottom_left_back_Node),
	//	bottom_right_front(bottom_right_front_Node),
	//	bottom_right_back(bottom_right_back_Node) {}
OctreeNode::OctreeNode() {};
OctreeNode::OctreeNode(double xminValue, double xmaxValue, double yminValue, double ymaxValue, double zminValue, double zmaxValue,int _tube_id)
{
	tube_id = _tube_id;
	xmin = xminValue;
	xmax = xmaxValue;
	ymin = yminValue;
	ymax = ymaxValue;
	zmin = zminValue;
	zmax = zmaxValue;
	/*top_left_front = top_left_back = top_right_front = top_right_back = bottom_left_front = bottom_left_back = bottom_right_front = bottom_right_back = NULL;*/
	children.resize(8);
	for (size_t i = 0; i < 8; ++i)
	{
		children[i] = NULL;
	}
}




//创建八叉树
void createOctree(OctreeNode* root,  double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
	//cout<<"处理中，请稍候……"<<endl;  
	//maxdepth = maxdepth - 1; //每递归一次就将最大递归深度-1  
	//if (maxdepth >= 0)
	//{
		root = new OctreeNode();
		//cout << "请输入节点值:";
		//root->data =9;//为节点赋值，可以存储节点信息，如物体可见性。由于是简单实现八叉树功能，简单赋值为9。  
		//cin >> root->data;  //为节点赋值  
		root->xmin = xmin; //为节点坐标赋值
		root->xmax = xmax;
		root->ymin = ymin;
		root->ymax = ymax;
		root->zmin = zmin;
		root->zmax = zmax;
		//double xm = (xmax + xmin) / 2;//计算节点个维度上的中间点
		//double ym = (ymax + ymin) / 2;
		//double zm = (zmax + zmin) / 2;
	//}
}
//int i = 1;
//先序遍历八叉树  
//template <class T>
//void preOrder(OctreeNode<T>*& p)
//{
//	if (p)
//	{
//		cout << i << ".当前节点的值为：" << p->data << "\n坐标为：";
//		cout << "xmin: " << p->xmin << " xmax: " << p->xmax;
//		cout << "ymin: " << p->ymin << " ymax: " << p->ymax;
//		cout << "zmin: " << p->zmin << " zmax: " << p->zmax;
//		i += 1;
//		cout << endl;
//		preOrder(p->top_left_front);
//		preOrder(p->top_left_back);
//		preOrder(p->top_right_front);
//		preOrder(p->top_right_back);
//		preOrder(p->bottom_left_front);
//		preOrder(p->bottom_left_back);
//		preOrder(p->bottom_right_front);
//		preOrder(p->bottom_right_back);
//		cout << endl;
//	}
//}
//求八叉树的深度  
//template<class T>
//int depth(OctreeNode<T>*& p)
//{
//	if (p == NULL)
//		return -1;
//	int h = depth(p->top_left_front);
//	return h + 1;
//}
//template<class T>
//int num(OctreeNode<T>*& p)
//{
//	if (p == NULL)
//		return 0;
//	return 1 + num(p->top_left_front) + num(p->top_left_back) + num(p->top_right_back) + num(p->top_right_front) + num(p->bottom_left_back) + num(p->bottom_left_front) + num(p->bottom_right_back) + num(p->bottom_right_front);
//}
//计算单位长度，为查找点做准备  
int cal(int num)
{
	int result = 1;
	if (1 == num)
		result = 1;
	else
	{
		for (int i = 1; i < num; i++)
			result = 2 * result;
	}
	return result;
}
//查找点  
//int maxdepth = 0;
//int times = 0;
//static double xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;
//int tmaxdepth = 0;
//double txm = 1, tym = 1, tzm = 1;
void find(OctreeNode* p, double x, double y, double z)
{
	double xm = (p->xmax - p->xmin) / 2;
	double ym = (p->ymax - p->ymin) / 2;
	double zm = (p->ymax - p->ymin) / 2;
	//times++;
	//if (x > xmax || x<xmin || y>ymax || y<ymin || z>zmax || z < zmin)
	//{
	//	/*cout << "该点不在场景中！" << endl;*/
	//	return;
	//}
	//if (x <= p->xmin + txm && x >= p->xmax - txm && y <= p->ymin + tym && y >= p->ymax - tym && z <= p->zmin + tzm && z >= p->zmax - tzm)
		//这里的判断主要是为了判断点是否处于一个更小划分区间里面
	//{
		//cout << endl << "找到该点！" << "该点位于" << endl;
		/*cout << "xmin: " << p->xmin << " xmax: " << p->xmax;
		cout << "ymin: " << p->ymin << " ymax: " << p->ymax;
		cout << "zmin: " << p->zmin << " zmax: " << p->zmax;
		cout << "节点内！" << endl;
		cout << "共经过" << times << "次递归！" << endl;*/
	//}
	if (x<(p->xmax - xm) && y>(p->ymax - ym) && z > (p->zmax - zm))
	{
		/*cout << "当前经过节点坐标：" << endl;
		cout << "xmin: " << p->xmin << " xmax: " << p->xmax;
		cout << "ymin: " << p->ymin << " ymax: " << p->ymax;
		cout << "zmin: " << p->zmin << " zmax: " << p->zmax;
		cout << endl;*/
		find(p->children[0], x, y, z);
	}
	else if (x < (p->xmax - xm) && y<(p->ymax - ym) && z>(p->zmax - zm))
	{
		//cout << "当前经过节点坐标：" << endl;
		/*cout << "xmin: " << p->xmin << " xmax: " << p->xmax;
		cout << "ymin: " << p->ymin << " ymax: " << p->ymax;
		cout << "zmin: " << p->zmin << " zmax: " << p->zmax;
		cout << endl;*/
		find(p->children[1], x, y, z);
	}
	else if (x > (p->xmax - xm) && y > (p->ymax - ym) && z > (p->zmax - zm))
	{
		/*cout << "当前经过节点坐标：" << endl;
		cout << "xmin: " << p->xmin << " xmax: " << p->xmax;
		cout << "ymin: " << p->ymin << " ymax: " << p->ymax;
		cout << "zmin: " << p->zmin << " zmax: " << p->zmax;
		cout << endl;*/
		find(p->children[2], x, y, z);
	}
	else if (x > (p->xmax - xm) && y<(p->ymax - ym) && z>(p->zmax - zm))
	{
		/*cout << "当前经过节点坐标：" << endl;
		cout << "xmin: " << p->xmin << " xmax: " << p->xmax;
		cout << "ymin: " << p->ymin << " ymax: " << p->ymax;
		cout << "zmin: " << p->zmin << " zmax: " << p->zmax;
		cout << endl;*/
		find(p->children[3], x, y, z);
	}
	else if (x<(p->xmax - xm) && y>(p->ymax - ym) && z < (p->zmax - zm))
	{
		/*cout << "当前经过节点坐标：" << endl;
		cout << "xmin: " << p->xmin << " xmax: " << p->xmax;
		cout << "ymin: " << p->ymin << " ymax: " << p->ymax;
		cout << "zmin: " << p->zmin << " zmax: " << p->zmax;
		cout << endl;*/
		find(p->children[4], x, y, z);
	}
	else if (x < (p->xmax - xm) && y < (p->ymax - ym) && z < (p->zmax - zm))
	{
		//cout << "当前经过节点坐标：" << endl;
		/*cout << "xmin: " << p->xmin << " xmax: " << p->xmax;
		cout << "ymin: " << p->ymin << " ymax: " << p->ymax;
		cout << "zmin: " << p->zmin << " zmax: " << p->zmax;
		cout << endl;*/
		find(p->children[5], x, y, z);
	}
	else if (x > (p->xmax - xm) && y > (p->ymax - ym) && z < (p->zmax - zm))
	{
		/*cout << "当前经过节点坐标：" << endl;
		cout << "xmin: " << p->xmin << " xmax: " << p->xmax;
		cout << "ymin: " << p->ymin << " ymax: " << p->ymax;
		cout << "zmin: " << p->zmin << " zmax: " << p->zmax;
		cout << endl;*/
		find(p->children[6], x, y, z);
	}
	else if (x > (p->xmax - xm) && y < (p->ymax - ym) && z < (p->zmax - zm))
	{
		/*cout << "当前经过节点坐标：" << endl;
		cout << "xmin: " << p->xmin << " xmax: " << p->xmax;
		cout << "ymin: " << p->ymin << " ymax: " << p->ymax;
		cout << "zmin: " << p->zmin << " zmax: " << p->zmax;
		cout << endl;*/
		find(p->children[7], x, y, z);
	}
}

bool canCover(OctreeNode* root, OctreeNode* tube_oct)
{
	if (root->xmin <= tube_oct->xmin && root->xmax >= tube_oct->xmax
		&& root->ymin <= tube_oct->ymin && root->ymax >= tube_oct->ymax
		&& root->zmin <= tube_oct->zmin && root->zmax >= tube_oct->zmax)
	{
		return true;
	}
	else {
		return false;
	}
}

void insert_tube(OctreeNode*root,  OctreeNode* tube_oct)
{
	double xmin = root->xmin;
	double xmax = root->xmax;
	double ymin = root->ymin;
	double ymax = root->ymax;
	double zmin = root->zmin;
	double zmax = root->zmax;

	double xm = (xmax +xmin) / 2;
	double ym = (ymax + ymin) / 2;
	double zm = (zmax + zmin) / 2;

	
	bool work = false;
	for (int i = 0; i < 8; ++i)
	{
		if (root->children[i] == NULL)
		{
			if (i == 0)
			{
				createOctree(root->children[i], xmin, xm, ym, ymax, zm, zmax); // top_left_front
			}
			if (i == 1)
			{
				createOctree(root->children[i], xmin, xm, ymin, ym, zm, zmax); // top_left_back
			}
			if (i == 2)
			{
				createOctree(root->children[i], xm, xmax, ym, ymax, zm, zmax); // top_right_front
			}
			if (i == 3)
			{
				createOctree(root->children[i], xm, xmax, ymin, ym, zm, zmax); // top_right_back
			}
			if (i == 4)
			{
				createOctree(root->children[i], xmin, xm, ym, ymax, zmin, zm); // bottom_left_front
			}
			if (i == 5)
			{
				createOctree(root->children[i], xmin, xm, ymin, ym, zmin, zm); // bottom_left_back
			}
			if (i == 6)
			{
				createOctree(root->children[i], xm, xmax, ym, ymax, zmin, zm); // bottom_right_front
			}
			if (i == 7)
			{
				createOctree(root->children[i], xm, xmax, ymin, ym, zmin, zm); // bottom_right_back
			}
		}

		if (canCover(root->children[i], tube_oct))
		{
			work = true;
			insert_tube(root->children[i], tube_oct);
			break;
		}
	}

	if (!work)
	{
		cv::Vec6d cuboid;
		cuboid[0] = tube_oct->xmin;
		cuboid[1] = tube_oct->xmax;
		cuboid[2] = tube_oct->ymin;
		cuboid[3] = tube_oct->ymax;
		cuboid[4] = tube_oct->zmin;
		cuboid[5] = tube_oct->zmax;

		root->tubes_cub.push_back(cuboid);
		root->tubes_oct.push_back(tube_oct);
	}
}

//int main()
//{
//	OctreeNode<double>* rootNode = NULL;
//	int choiced = 0;
//	cout << "系统开始前请先创建八叉树" << endl;
//	cout << "请输入最大递归深度：" << endl;
//	cin >> maxdepth;
//	cout << "请输入外包盒坐标，顺序如下：xmin,xmax,ymin,ymax,zmin,zmax" << endl;
//	cin >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
//	if (maxdepth >= 0 || xmax > xmin || ymax > ymin || zmax > zmin || xmin > 0 || ymin > 0 || zmin > 0)
//	{
//		tmaxdepth = cal(maxdepth);
//		txm = (xmax - xmin) / tmaxdepth;
//		tym = (ymax - ymin) / tmaxdepth;
//		tzm = (zmax - zmin) / tmaxdepth;
//		createOctree(rootNode, maxdepth, xmin, xmax, ymin, ymax, zmin, zmax);
//	}
//	while (true)
//	{
//		system("cls");
//
//		cout << "请选择操作：\n";
//		cout << "\t1.计算空间中区域的个数\n";
//		cout << "\t2.先序遍历八叉树\n";
//		cout << "\t3.查看树深度\n";
//		cout << "\t4.查找节点   \n";
//		cout << "\t0.退出\n";
//		cin >> choiced;
//
//		if (choiced == 0)
//			return 0;
//		if (choiced == 1)
//		{
//			system("cls");
//			cout << "空间区域个数" << endl;
//			cout << num(rootNode);
//		}
//
//		if (choiced == 2)
//		{
//			system("cls");
//			cout << "先序遍历八叉树结果：/n";
//			i = 1;
//			preOrder(rootNode);
//			cout << endl;
//			system("pause");
//		}
//		if (choiced == 3)
//		{
//			system("cls");
//			int dep = depth(rootNode);
//			cout << "此八叉树的深度为" << dep + 1 << endl;
//			system("pause");
//		}
//		if (choiced == 4)
//		{
//			system("cls");
//			cout << "请输入您希望查找的点的坐标，顺序如下：x,y,z\n";
//			double x, y, z;
//			cin >> x >> y >> z;
//			times = 0;
//			cout << endl << "开始搜寻该点……" << endl;
//			find(rootNode, x, y, z);
//			system("pause");
//		}
//		else
//		{
//			system("cls");
//			cout << "\n\n错误选择！\n";
//			system("pause");
//		}
//	}


bool cmp(pair<int, cv::Vec6d>& a, pair<int, cv::Vec6d>& b)
{
	return a.second[5] - a.second[4] > b.second[5] - b.second[4];
}

bool bb_overlap(cv::Vec6d& a, cv::Vec6d& b)
{
	if (a[0] >= b[1])
	{
		return false;
	}
		
	if (a[2] >= b[3])
	{
		return false;
	}
	
	if (a[1] <= b[0]) 
	{
		return false;
	}

	if (a[3] <= b[2]) 
	{
		return false;
	}

	return true;
}

//double cal_colli_area(cv::Vec6d& a, cv::Vec6d& b)
//{
//	/*if (a[0] != b[0])
//	{
//		return 0.0;
//	}*/
//	double x1 = a[1], y1 = a[2], w1 = a[3], h1 = a[4];
//	double x2 = b[1], y2 = b[2], w2 = b[3], h2 = b[4];
//
//	double endx = max(x1 + w1, x2 + w2);
//	double startx = min(x1, x2);
//	double width = w1 + w2 - (endx - startx);  // 重叠部分宽
//	double endy = max(y1 + h1, y2 + h2);
//	double starty = min(y1, y2);
//	double height = h1 + h2 - (endy - starty);  // 重叠部分高
//
//
//	double area1 = w1 * h1;
//	double area2 = w2 * h2;
//
//	if (width > 0 && height > 0) {
//		//int area = width * height;  // 重叠部分面积
//		return width * height;
//
//	}
//	else {
//		// 不重叠：算出来的width或height小于等于0
//		return 0.0;
//	}
//}

int insert(std::map<int, cv::Vec6d>& syn_tube_cuboid_vec, pair<int, cv::Vec6d> insert_tube, Context &context, int &beta)
{
	int res_frame_idx = -1;
	int minColli = INT_MAX;
	int tube_len = insert_tube.second[5] - insert_tube.second[4] + 1.0;
	// 在所有可行的时间位置上判断
	for (int beg_frame_idx = 0; ;++beg_frame_idx) // beg_frame_idx <= context.synoplen - tube_len
	{
		int colli = 0;

		int end_frame_idx = beg_frame_idx + tube_len - 1;
		for (auto it :syn_tube_cuboid_vec)
		{
			int syn_tube_beg_frame_idx = it.second[4];
			int syn_tube_end_frame_idx = it.second[5];

			int shared_frame_beg_idx = std::max(beg_frame_idx, syn_tube_beg_frame_idx);
			int shared_frame_end_idx = std::min(end_frame_idx, syn_tube_end_frame_idx);

			/*for (int shared_frame_idx = shared_frame_beg_idx; shared_frame_idx <= shared_frame_end_idx; ++shared_frame_idx)
			{

			}*/
			// 时间和空间上都有重叠
			/*if (shared_frame_end_idx >= shared_frame_beg_idx && bb_overlap(tube_cuboid, syn_tube_cuboid_vec[j].second))
			{
				colli += shared_frame_end_idx - shared_frame_beg_idx + 1;
			}*/
			for (int shared_frame_idx = shared_frame_beg_idx; shared_frame_idx <= shared_frame_end_idx; ++shared_frame_idx)
			{
				cv::Vec6d bbox1 = context.bbox[it.first-1][shared_frame_idx - syn_tube_beg_frame_idx];
				cv::Vec6d cuboid1;
				cuboid1[0] = bbox1[1];
				cuboid1[1] = bbox1[1] + bbox1[3];
				cuboid1[2] = bbox1[2];
				cuboid1[3] = bbox1[2] + bbox1[4];
				cuboid1[4] = cuboid1[5] = 0.0;

				cv::Vec6d bbox2 = context.bbox[insert_tube.first-1][shared_frame_idx - beg_frame_idx];
				cv::Vec6d cuboid2;
				cuboid2[0] = bbox2[1];
				cuboid2[1] = bbox2[1] + bbox2[3];
				cuboid2[2] = bbox2[2];
				cuboid2[3] = bbox2[2] + bbox2[4];
				cuboid2[4] = cuboid2[5] = 0.0;
				// 判断是否碰撞
				if (bb_overlap(cuboid1, cuboid2))
				{
					colli += 1;
				}
			}
		}

		// 找到一个具有更少碰撞的时间位置
		if (colli < beta)
		{
			res_frame_idx = beg_frame_idx;
			break;
		}
	}

	// 插入tube
	cv::Vec6d tube_res_cuboid = insert_tube.second;

	tube_res_cuboid[4] = double(res_frame_idx);
	tube_res_cuboid[5] = double(res_frame_idx + tube_len - 1);

	syn_tube_cuboid_vec[insert_tube.first] = tube_res_cuboid;
	return res_frame_idx;
}


std::map<int, cv::Vec6d> refine(std::map<int, cv::Vec6d>syn_tube_cuboid_vec, pair<int, cv::Vec6d> insert_tube, Context &context, int &beta, int &syn_last_frame_idx)
{
	// 先把insert_tube 移除
	syn_tube_cuboid_vec.erase(insert_tube.first);
	// 找到所有需要重新insert的tube
	std::vector<pair<int, cv::Vec6d>>refine_tubes_cuboid;

	auto it = syn_tube_cuboid_vec.begin();
	while (it != syn_tube_cuboid_vec.end())
	{
		int syn_tube_id = it->first;
		cv::Vec6d syn_tube_cuboid = it->second;
		if (syn_tube_cuboid[5]-syn_tube_cuboid[4] < insert_tube.second[5]-insert_tube.second[4] && bb_overlap(it->second, insert_tube.second))
		{
			refine_tubes_cuboid.push_back(pair<int, cv::Vec6d>(syn_tube_id, syn_tube_cuboid));
			// 先把与该insert_tube碰撞的tube移除
			syn_tube_cuboid_vec.erase(it++);
		}
		else {
			++it;
		}
	}
	
	// 插入insert_tube
	int res_frame_idx=insert(syn_tube_cuboid_vec, insert_tube, context, beta);
	
	int tmp_syn_last_frame_idx = max(res_frame_idx + int(insert_tube.second[5] - insert_tube.second[4]), syn_last_frame_idx);

	if (tmp_syn_last_frame_idx > syn_last_frame_idx)
	{
		std::map<int, cv::Vec6d>tmp;
		tmp.clear();
		return tmp;
	}

	// 对需要refine的tube 按 tube_len 降序排序
	sort(refine_tubes_cuboid.begin(), refine_tubes_cuboid.end(), cmp);

	for (size_t i = 0; i < refine_tubes_cuboid.size(); ++i)
	{
		int res_frame_idx_for_refine = insert(syn_tube_cuboid_vec, refine_tubes_cuboid[i], context, beta);
		int tube_len = refine_tubes_cuboid[i].second[5] - refine_tubes_cuboid[i].second[4] + 1.0;

		if (res_frame_idx_for_refine + tube_len - 1 > syn_last_frame_idx)
		{
			std::map<int, cv::Vec6d>tmp;
			tmp.clear();
			return tmp;
		}
		tmp_syn_last_frame_idx = max(tmp_syn_last_frame_idx, res_frame_idx_for_refine + tube_len - 1);
	}
	
	syn_last_frame_idx = tmp_syn_last_frame_idx;
	return syn_tube_cuboid_vec;
}

void octree_tube_rerrange(Context& context, std::map<int, int>& tubes_group_id, double speed)
{
	// 碰撞阈值
	int beta = 400;
	
	// 摘要长度阈值
	//int synoplen_threshold = 1000;

	int syn_last_frame_idx = 0;

	context.tube_num = 193 - 1 + 1;
	context.scale = 1.0;
    // 获取 bbox and occlu 信息
	read_boundingbox(context.filepath + "/boundingbox", context.bbox, context.tube_num, context.scale, context.video_height, context.video_width);
	//read_occlu(context.filepath, context.orig_colli);

	std::cout << "read complete." << std::endl;
	// 加速
	for (int i = 0; i < int(context.bbox.size()); ++i)
	{
		int beg = context.bbox[i][0][0];
		int tmp = int(context.bbox[i].size()) - 1;
		int end = context.bbox[i][tmp][0];
		// 计算加速后的新结束帧号
		int newEnd = max(beg, beg + int(round(speed * (end - beg + 1))) - 1);
		std::vector<cv::Vec6d> bboxes_cur_tube;
		// 获取加速后，当前tube的各帧的bbox
		for (int kk = beg; kk <= newEnd; ++kk)
		{
			int k;
			if (beg == newEnd)
			{
				k = beg;
			}
			else
			{
				k = (int)(((double)(kk - beg) / (newEnd - beg)) * (end - beg) + beg);
			}
			// 更新帧号
			cv::Vec6d bbox = context.bbox[i][k - beg];
			bbox[0] = kk;
			// 记录原始帧号
			bbox[5] = k;
			bboxes_cur_tube.push_back(bbox);
		}
		context.bbox[i].assign(bboxes_cur_tube.begin(), bboxes_cur_tube.end());
	}


	// 构建全局oct
	//OctreeNode* global_tree = new OctreeNode(0,context.video_width,0,context.video_height,0,context.synoplen);

	std::vector<pair<int,cv::Vec6d>> tube_cuboid_vec;

	int tube_num_longer_synoplen = 0;
	for (size_t i = 0; i < context.bbox.size(); ++i)
	{
		double xmin = 1e5, xmax = 0.0;
		double ymin = 1e5, ymax = 0.0;

		int tmp = int(context.bbox[i].size()) - 1;
		double zmin = context.bbox[i][0][0], zmax = context.bbox[i][tmp][0];


		for (size_t j = 0; j < context.bbox[i].size(); ++j)
		{
			xmin = min(xmin, context.bbox[i][j][1]);
			xmax = max(xmax, context.bbox[i][j][1] + context.bbox[i][j][3]);

			ymin = min(ymin, context.bbox[i][j][2]);
			ymax = max(ymax, context.bbox[i][j][2] + context.bbox[i][j][4]);
		}

		// 构建tube_oct，并insert
		/*OctreeNode* tube_tree = new OctreeNode(xmin, xmax, ymin, ymax, zmin, zmax, i + 1);
		insert_tube(global_tree, tube_tree);*/
		cv::Vec6d cuboid;
		cuboid[0] = xmin;
		cuboid[1] = xmax;
		cuboid[2] = ymin;
		cuboid[3] = ymax;
		cuboid[4] = zmin;
		cuboid[5] = zmax;

		pair<int, cv::Vec6d>cuboid_with_id;
		cuboid_with_id.first = i + 1;
		cuboid_with_id.second = cuboid;
		tube_cuboid_vec.push_back(cuboid_with_id);

		if (int(zmax - zmin + 1.0) > context.synoplen)
		{
			++tube_num_longer_synoplen;
		}
	}

	std::cout << "tube_num_longer_synoplen: " << tube_num_longer_synoplen << std::endl;

	for (size_t i = 0; i < tube_cuboid_vec.size(); ++i)
	{
		syn_last_frame_idx = tube_cuboid_vec[i].second[5] - tube_cuboid_vec[i].second[4];
		
		if(syn_last_frame_idx <= context.synoplen-1)
		{
			break;
		}
	}
	
	// 摘要空间中的已经插入的tube cuboid
	std::map<int, cv::Vec6d>syn_tube_cuboid_vec;
	
	std::vector<pair<int,cv::Vec6d>> tmp_tube_cuboid_vec;
	tmp_tube_cuboid_vec.assign(tube_cuboid_vec.begin(), tube_cuboid_vec.end());

	while (true)
	{
		// 目前无法插入的tube
		std::vector<pair<int,cv::Vec6d>> tube_dropped;

		for (size_t i = 0; i < tube_cuboid_vec.size(); ++i)
		{
			int tube_id = tube_cuboid_vec[i].first;
			cv::Vec6d tube_cuboid = tube_cuboid_vec[i].second;

			int tube_len = tube_cuboid[5] - tube_cuboid[4] + 1.0;
			if (tube_len > context.synoplen)
			{
				continue;
			}

			// 插入
			int res_frame_idx = insert(syn_tube_cuboid_vec, pair<int, cv::Vec6d>(tube_id, tube_cuboid), context, beta);

			// 计算当前tube的结束位置
			int end_frame_idx = res_frame_idx + tube_len - 1;

			// refine
			if (end_frame_idx > syn_last_frame_idx)
			{
				std::map<int, cv::Vec6d>refine_res;
				refine_res = refine(syn_tube_cuboid_vec, pair<int, cv::Vec6d>(tube_id, tube_cuboid), context, beta, syn_last_frame_idx);
				// refine success
				if (!refine_res.empty())
				{
					syn_tube_cuboid_vec = refine_res;
				}
				// refine fail
				else
				{
					// 根据摘要长度阈值限制tube数量
					if (end_frame_idx + 1 > context.synoplen)
					{
						syn_tube_cuboid_vec.erase(tube_id);

						tube_dropped.push_back(tube_cuboid_vec[i]);
						continue;
					}

					syn_last_frame_idx = max(syn_last_frame_idx, end_frame_idx);
				}
			}
			// dont need to refine
			else
			{
				// 根据摘要长度阈值限制tube数量
				if (end_frame_idx + 1 > context.synoplen)
				{
					syn_tube_cuboid_vec.erase(tube_id);

					tube_dropped.push_back(tube_cuboid_vec[i]);
					continue;
				}

				syn_last_frame_idx = max(syn_last_frame_idx, end_frame_idx);
			}

			std::cout << "tube " << tube_id << " complete!" << std::endl;
		}

		if (tube_dropped.empty() || tube_cuboid_vec.size() == tube_dropped.size())
		{
			break;
		}

		tube_cuboid_vec.assign(tube_dropped.begin(), tube_dropped.end());
		std::cout << "one step rearrange complete" << std::endl;
	}
	

	// 更新摘要长度
	//context.synoplen = syn_last_frame_idx + 1;
	//---------------------------------------------------------
	//context.synoplen = synoplen_threshold;

	std::cout << syn_tube_cuboid_vec.size() << std::endl;
	// 写入 all_sges和bestCopy
	std::map<int, std::vector<Segment*>> best;
	int max_bb = 0;
	for (auto it: syn_tube_cuboid_vec)
	{
		int tube_id = it.first;
		int a = int(tmp_tube_cuboid_vec[tube_id - 1].second[4]);
		int b = int(tmp_tube_cuboid_vec[tube_id - 1].second[5]);

		std::vector<Segment*> tmp_seg;
		tmp_seg.push_back(new Segment(a, b, int(it.second[4]), int(it.second[5]), tube_id - 1, b - a + 1));

		max_bb = max(max_bb, int(it.second[5]));

		best[tube_id] = tmp_seg;
	}

	std::cout << "max_bb: " << max_bb << std::endl;

	for (int tube_id = 1; tube_id <= context.tube_num; ++tube_id)
	{
		if (best.count(tube_id) > 0)
		{
			context.all_segs.push_back(best[tube_id]);
		}
		else 
		{
			std::vector<Segment*> tmp_seg;
			tmp_seg.clear();
			context.all_segs.push_back(tmp_seg);
		}
	}

	for (size_t i = 1; i <= context.tube_num; ++i)
	{
		tubes_group_id[i] = i;
	}

	// 计算碰撞
	/*std::map<int,std::vector<cv::Vec6d>>bbox_per_frame;
	for (size_t i = 0; i < context.all_segs.size(); ++i)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); ++j)
		{
			Segment* seg = context.all_segs[i][j];
			for (int frame_idx = seg->aa_; frame_idx <= seg->bb_; ++frame_idx)
			{
				cv::Vec6d patch = context.bbox[i][frame_idx - seg->aa_];
				bbox_per_frame[frame_idx].push_back(patch);
			}
		}
	}

	double colli_area = 0.0;
	for (auto it : bbox_per_frame)
	{
		for (int i = 0; i < int(it.second.size())-1; ++i)
		{
			for (int j = i + 1; j < int(it.second.size()); ++j)
			{
				colli_area += cal_colli_area(it.second[i], it.second[j]);
			}
		}
	}*/

	std::ofstream ofile(context.filepath + "/metrics.txt", ofstream::app);
	//ofile << "colli_area: "<< colli_area << "\n";
	ofile << "beta: " << beta << "\n";
	ofile << "number of dropped tube: " << context.tube_num - int(syn_tube_cuboid_vec.size()) << "\n";

	std::cout << "number of dropped tube: " << context.tube_num - int(syn_tube_cuboid_vec.size()) << std::endl;
}