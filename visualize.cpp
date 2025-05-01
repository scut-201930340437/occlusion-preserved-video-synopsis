#include "visualize.h"
#include "frame_seq_fusion.h"
#include "energy.h"
#include <random>
#include <iostream>
#include <fstream>
#include <cmath>

/* Use graph::dfs_connected_graph_ofA instead */
//void GetAllSegmentsFromRoot(Segment* root, std::vector<Segment*>& segs)
//{
//	segs.push_back(root);
//	root->is_traversed = true;
//
//	EdgeNode* edge = root->firstedge;
//
//	while (edge != NULL)
//	{
//		if (edge->is_retained && edge->seg->is_traversed == false)
//			GetAllSegmentsFromRoot(edge->seg, segs);
//		edge = edge->next;
//	}
//}


cv::Rect get_bbox_of_obj_by_frame(std::vector<std::vector<cv::Vec6d>>& bbox, int tubeid, int frameid)
{
	/*for (size_t i = 0; i < bbox[tubeid].size(); i++)
	{
		if (bbox[tubeid][i][0] == frameid)
		{
			return cv::Rect(int(bbox[tubeid][i][1]), int(bbox[tubeid][i][2]), int(bbox[tubeid][i][3]), int(bbox[tubeid][i][4]));
		}
	}*/
	int beg = bbox[tubeid][0][0];
	int tmp = bbox[tubeid].size() - 1;
	int end = bbox[tubeid][tmp][0];
	if (frameid < beg || frameid > end)
	{
		std::cout << "tube: " << tubeid << std::endl;
		std::cout << "frame_id: " << frameid << std::endl;
		std::cout << "get_rect_of_obj_by_frame throws exception" << std::endl;
		return cv::Rect();
	}
	else 
	{
		int i = frameid - beg;
		return cv::Rect(int(bbox[tubeid][i][1]), int(bbox[tubeid][i][2]), int(bbox[tubeid][i][3]), int(bbox[tubeid][i][4]));
	}
	
}

struct VisMeta
{
	int src_frame_id;
	int dst_frame_id;
	int object_id;
	double scale;
	cv::Rect box; 
};


cv::Scalar random_color(int id)
{
	std::default_random_engine e;
	std::uniform_int_distribution<int> u(0, 255); // 左闭右闭区间
	e.seed(id);
	int r = u(e);
	int g = u(e);
	int b = u(e);

	/*srand(group_id*5);
	int r = rand() % 256;
	int g = rand() % 256;
	int b = rand() % 256;*/

	return cv::Scalar(r, g, b);
}

void draw_dotted_line2(cv::Mat &img, cv::Point2d p1, cv::Point2d p2, cv::Scalar color, int thickness=1)
{
	double n = 5; //线长度
	double w = p2.x - p1.x, h = p2.y - p1.y;
	double l = sqrtf(w * w + h * h);
	// 矫正线长度，使线个数为奇数
	int m = l / n;
	m = m % 2 ? m : m + 1;
	n = l / m;

	circle(img, p1, 1, color, thickness); // 画起点
	circle(img, p2, 1, color, thickness); // 画终点
	// 画中间点
	if (p1.y == p2.y) //水平线：y = m
	{
		double x1 = min(p1.x, p2.x);
		double x2 = max(p1.x, p2.x);
		for (double x = x1, n1 = 2 * n; x < x2; x = x + n1)
			line(img, cv::Point2d(x, p1.y), cv::Point2d(x + n, p1.y), color, thickness);
	}
	else if (p1.x == p2.x) //垂直线, x = m
	{
		double y1 = min(p1.y, p2.y);
		double y2 = max(p1.y, p2.y);
		for (double y = y1, n1 = 2 * n; y < y2; y = y + n1)
			line(img, cv::Point2d(p1.x, y), cv::Point2d(p1.x, y + n), color, thickness);
	}
	else // 倾斜线，与x轴、y轴都不垂直或平行
	{
		// 直线方程的两点式：(y-y1)/(y2-y1)=(x-x1)/(x2-x1) -> y = (y2-y1)*(x-x1)/(x2-x1)+y1
		double n1 = n * abs(w) / l;
		double k = h / w;
		double x1 = min(p1.x, p2.x);
		double x2 = max(p1.x, p2.x);
		for (double x = x1, n2 = 2 * n1; x < x2; x = x + n2)
		{
			cv::Point2d p3 = cv::Point2d(x, k * (x - p1.x) + p1.y);
			cv::Point2d p4 = cv::Point2d(x + n1, k * (x + n1 - p1.x) + p1.y);
			line(img, p3, p4, color, thickness);
		}
	}
}


void visualize2(Context& context, std::vector<std::vector<Segment*>> best)
{

	int synlen = context.synoplen;

	cv::VideoCapture cap(context.filepath + "/src.mp4");


	int src_video_len = int(cap.get(cv::CAP_PROP_FRAME_COUNT));
	double fps = cap.get(cv::CAP_PROP_FPS);

	std::vector<cv::Mat> frames;
	cv::Mat tmp2 = cv::imread(context.filepath + "/background.png");
	cv::resize(tmp2, tmp2, cv::Size(), context.scale, context.scale);
	int video_height = tmp2.rows;
	int video_width = tmp2.cols;

	cv::VideoWriter capw;

	if (!capw.open(context.filepath + "/synopsis_step2.mp4", cv::VideoWriter::fourcc('H', '2', '6', '4'), fps, tmp2.size()))
	{
		std::cout << "Can not write!" << std::endl;
		return;
	}

	for (int i = 0; i < synlen /*context.synoplen*/; i++)
	{
		frames.push_back(tmp2.clone());
	}

	std::vector<std::vector<VisMeta>> metas;
	metas.resize(src_video_len);

	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		/*if (context.circulated_tube[i] == 1)
		{
			continue;
		}*/

		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			//std::cout << i << "," << j << std::endl;

			int aa = context.all_segs[i][j]->aa_;
			int bb = context.all_segs[i][j]->bb_;
			double scale = context.all_segs[i][j]->scale_;


			for (int kk = aa; kk <= bb; kk++)
			{
				if (kk >= 0 && kk < synlen)
				{
					int k = context.all_segs[i][j]->interp_frame_No(kk);
					cv::Rect roi = get_bbox_of_obj_by_frame(context.bbox, i, k);


					VisMeta vm;
					vm.src_frame_id = k; //原视频中的帧
					vm.dst_frame_id = kk; //摘要视频中的帧
					vm.scale = scale;
					vm.box = roi;
					vm.object_id = i + 1;
					metas[k].push_back(vm);


				}
			}
		}
	}
	//cap.set(cv::CAP_PROP_POS_FRAMES, 0);


	for (int i = 0; i < src_video_len; i++) //按照原视频中的帧
	{
		cv::Mat iimg;
		cap.read(iimg);

		if (i % 1000 == 0)
		{
			std::cout << i << "/" << src_video_len << std::endl;
		}

		//std::cout << "a" << std::endl;


		if (metas[i].empty())
		{
			continue;
		}
		//std::cout << "b" << std::endl;

		cv::resize(iimg, iimg, cv::Size(), context.scale, context.scale);
		//std::cout << "c" << std::endl;

		//std::cout << iimg.rows << iimg.cols << std::endl;

		for (size_t j = 0; j < metas[i].size(); j++)
		{
			VisMeta vm = metas[i][j];
			cv::Rect roi = vm.box;
			int dx = int(roi.width * (vm.scale - 1) / 2.0);
			int dy = int(roi.height * (vm.scale - 1) / 2.0);

			cv::Rect new_roi(roi.x - dx, roi.y - dy, roi.width + 2 * dx, roi.height + 2 * dy);

			if (roi.width <= 0 || roi.height <= 0 || roi.width > video_width || roi.height > video_height)
				continue;

			if (new_roi.width <= 0 || new_roi.height <= 0 || new_roi.width > video_width || new_roi.height > video_height)
				continue;

			//std::cout << roi << std::endl;

			cv::Mat mask(new_roi.height, new_roi.width, CV_8U, cv::Scalar(1));

			cv::Mat new_img;
			cv::resize(iimg(roi), new_img, cv::Size(new_roi.width, new_roi.height));//iimg为原视频中的帧



			new_img.copyTo(frames[vm.dst_frame_id](new_roi), mask); //dst摘要视频



		}

		for (size_t j = 0; j < metas[i].size(); j++)
		{
			VisMeta vm = metas[i][j];
			cv::Rect roi = vm.box;
			int dx = int(roi.width * (vm.scale - 1) / 2.0);
			int dy = int(roi.height * (vm.scale - 1) / 2.0);

			cv::Rect new_roi(roi.x - dx, roi.y - dy, roi.width + 2 * dx, roi.height + 2 * dy);

			if (roi.width <= 0 || roi.height <= 0 || roi.width > video_width || roi.height > video_height)
				continue;

			if (new_roi.width <= 0 || new_roi.height <= 0 || new_roi.width > video_width || new_roi.height > video_height)
				continue;

			//std::cout << roi << std::endl;

			cv::Mat mask(new_roi.height, new_roi.width, CV_8U, cv::Scalar(1));

			cv::Mat new_img;
			cv::resize(iimg(roi), new_img, cv::Size(new_roi.width, new_roi.height));


			//new_img.copyTo(frames[vm.dst_frame_id](new_roi), mask);

			//ID
			cv::rectangle(frames[vm.dst_frame_id], roi, cv::Scalar(0, 0, 255));
			char text[1024];
			sprintf_s(text, 1024, "%d", vm.object_id);
			cv::putText(frames[vm.dst_frame_id], text, cv::Point(roi.x /*+ roi.width / 4*/, roi.y - 5),
				cv::FONT_HERSHEY_PLAIN, 1, cv::Scalar(0, 0, 255));

		}

	}


	for (size_t i = 0; i < frames.size(); i++)
	{
		if (i % 1000 == 0)
			std::cout << i << std::endl;
		capw.write(frames[i]);
	}
	//std::cout << "frame number:"<< frames.size() << std::endl;
}

void visualize1(Context& context, std::vector<std::vector<Segment*>> best)
{

	int synlen = context.synoplen;

	cv::VideoCapture cap(context.filepath + "/src.mp4");


	int src_video_len = int(cap.get(cv::CAP_PROP_FRAME_COUNT));
	double fps = cap.get(cv::CAP_PROP_FPS);

	std::vector<cv::Mat> frames;
	cv::Mat tmp2 = cv::imread(context.filepath + "/background.png");
	cv::resize(tmp2, tmp2, cv::Size(), context.scale, context.scale);
	int video_height = tmp2.rows;
	int video_width = tmp2.cols;

	cv::VideoWriter capw;

	if (!capw.open(context.filepath + "/synopsis_step1.mp4", cv::VideoWriter::fourcc('H', '2', '6', '4'), fps, tmp2.size()))
	{
		std::cout << "Can not write!" << std::endl;
		return;
	}

	for (int i = 0; i < synlen /*context.synoplen*/; i++)
	{
		frames.push_back(tmp2.clone());
	}

	std::vector<std::vector<VisMeta>> metas;
	metas.resize(src_video_len);

	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		/*if (context.circulated_tube[i] == 1)
		{
			continue;
		}*/

		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			//std::cout << i << "," << j << std::endl;

			int aa = context.all_segs[i][j]->aa_;
			int bb = context.all_segs[i][j]->bb_;
			double scale = context.all_segs[i][j]->scale_;


			for (int kk = aa; kk <= bb; kk++)
			{
				if (kk >= 0 && kk < synlen)
				{
					int k = context.all_segs[i][j]->interp_frame_No(kk);
					cv::Rect roi = get_bbox_of_obj_by_frame(context.bbox, i, k);
					

					VisMeta vm;
					vm.src_frame_id = k; //原视频中的帧
					vm.dst_frame_id = kk; //摘要视频中的帧
					vm.scale = scale;
					vm.box = roi;
					vm.object_id = i + 1;
					metas[k].push_back(vm);

					
				}
			}
		}
	}
	//cap.set(cv::CAP_PROP_POS_FRAMES, 0);


	for (int i = 0; i < src_video_len; i++) //按照原视频中的帧
	{
		cv::Mat iimg;
		cap.read(iimg);

		if (i % 1000 == 0)
		{
			std::cout << i << "/" << src_video_len << std::endl;
		}

		//std::cout << "a" << std::endl;


		if (metas[i].empty())
		{
			continue;
		}
		//std::cout << "b" << std::endl;

		cv::resize(iimg, iimg, cv::Size(), context.scale, context.scale);
		//std::cout << "c" << std::endl;

		//std::cout << iimg.rows << iimg.cols << std::endl;

		for (size_t j = 0; j < metas[i].size(); j++)
		{
			VisMeta vm = metas[i][j];
			cv::Rect roi = vm.box;
			int dx = int(roi.width * (vm.scale - 1) / 2.0);
			int dy = int(roi.height * (vm.scale - 1) / 2.0);

			cv::Rect new_roi(roi.x - dx, roi.y - dy, roi.width + 2 * dx, roi.height + 2 * dy);

			if (roi.width <= 0 || roi.height <= 0 || roi.width > video_width || roi.height > video_height)
				continue;

			if (new_roi.width <= 0 || new_roi.height <= 0 || new_roi.width > video_width || new_roi.height > video_height)
				continue;

			//std::cout << roi << std::endl;

			cv::Mat mask(new_roi.height, new_roi.width, CV_8U, cv::Scalar(1));

			cv::Mat new_img;
			cv::resize(iimg(roi), new_img, cv::Size(new_roi.width, new_roi.height));//iimg为原视频中的帧

			

			new_img.copyTo(frames[vm.dst_frame_id](new_roi), mask); //dst摘要视频

			

		}

		for (size_t j = 0; j < metas[i].size(); j++)
		{
			VisMeta vm = metas[i][j];
			cv::Rect roi = vm.box;
			int dx = int(roi.width * (vm.scale - 1) / 2.0);
			int dy = int(roi.height * (vm.scale - 1) / 2.0);

			cv::Rect new_roi(roi.x - dx, roi.y - dy, roi.width + 2 * dx, roi.height + 2 * dy);

			if (roi.width <= 0 || roi.height <= 0 || roi.width > video_width || roi.height > video_height)
				continue;

			if (new_roi.width <= 0 || new_roi.height <= 0 || new_roi.width > video_width || new_roi.height > video_height)
				continue;

			//std::cout << roi << std::endl;

			cv::Mat mask(new_roi.height, new_roi.width, CV_8U, cv::Scalar(1));

			cv::Mat new_img;
			cv::resize(iimg(roi), new_img, cv::Size(new_roi.width, new_roi.height));


			//new_img.copyTo(frames[vm.dst_frame_id](new_roi), mask);

			//ID
			cv::rectangle(frames[vm.dst_frame_id], roi, cv::Scalar(0, 0, 255));
			char text[1024];
			sprintf_s(text, 1024, "%d", vm.object_id);
			cv::putText(frames[vm.dst_frame_id], text, cv::Point(roi.x /*+ roi.width / 4*/, roi.y - 5),
				cv::FONT_HERSHEY_PLAIN, 1, cv::Scalar(0, 0, 255));

		}

	}


	for (size_t i = 0; i < frames.size(); i++)
	{
		if (i % 1000 == 0)
			std::cout << i << std::endl;
		capw.write(frames[i]);
	}
	//std::cout << "frame number:"<< frames.size() << std::endl;
}

void visualize(Context& context, std::map<int, int>& tubes_group_id, double speed, int chrono_slack)
{
	int synlen = context.synoplen;

	cv::VideoCapture cap(context.filepath + "/src.mp4");
	int src_video_len = int(cap.get(cv::CAP_PROP_FRAME_COUNT));
	double fps = cap.get(cv::CAP_PROP_FPS);

	cv::Mat tmp2 = cv::imread(context.filepath + "/background.png");
	cv::resize(tmp2, tmp2, cv::Size(), context.scale, context.scale);
	int video_height = tmp2.rows;
	int video_width = tmp2.cols;

	cv::VideoWriter capw;
	
	if (!capw.open(context.filepath + "/synopsis.mp4", cv::VideoWriter::fourcc('m', 'p', 'v', '4'), fps, tmp2.size()))
	{
		std::cout << "Can not write!" << std::endl;
		return;
	}
	
	std::vector<std::vector<VisMeta>> metas;
	metas.resize(src_video_len);

	// 可视化的tube数量
	int visualize_tube_num = 0;
	context.visualize_state.resize(int(context.all_segs.size()), true);


	//std::vector<int>vis_tube_vec = { 2, 4, 8, 10,11,14,15,17,18,19,20 };
	/*for (size_t i = 0; i < vis_tube_vec.size(); ++i)
	{
		context.visualize_state[vis_tube_vec[i] - 1] = true;
	}*/

	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		/*if (context.circulated_tube[i] == 1)
		{
			continue;
		}*/


		/*if (not_vis_numS.count(i+tube_beg) > 0)
		{
			continue;
		}*/

		// 不可视化边缘出现 边缘消失的tube
		/*double xmin = 1000000000000.0, xmax = -1.0, ymin = 100000000000000.0, ymax = -1.0;
		for (size_t ij = 0; ij < context.bbox[i].size(); ++ij)
		{
			xmin = min(xmin, context.bbox[i][ij][1]);
			xmax = max(xmax, context.bbox[i][ij][1]);

			ymin = min(ymin, context.bbox[i][ij][2]);
			ymax = max(ymax, context.bbox[i][ij][2]);
		}

		if ((xmax - xmin) < 140.0 && (ymax - ymin) < 140.0)
		{
			continue;
		}*/

		int tmp = int(context.bbox[i].size())-1;
		
		/*if (context.bbox[i][0][2] > 290.0)
		{
			if (context.bbox[i][0][1] < 770.0 && context.bbox[i][0][1] > 160.0 && context.bbox[i][0][2] < 530.0)
			{
				continue;
			}
		}
		else
		{
			if (context.bbox[i][0][1] < 750.0 && context.bbox[i][0][1] > 160.0 && context.bbox[i][0][2] > 170.0)
			{
				continue;
			}
		}

		if (context.bbox[i][tmp][2] > 290.0)
		{
			if (context.bbox[i][tmp][1] < 770.0 && context.bbox[i][tmp][1] > 160.0 && context.bbox[i][tmp][2] < 530.0)
			{
				continue;
			}
		}
		else
		{
			if (context.bbox[i][tmp][1] < 750.0 && context.bbox[i][tmp][1] > 160.0 && context.bbox[i][tmp][2] > 170.0)
			{
				continue;
			}
		}*/


		/*if (context.bbox[i][0][2] > 350.0)
		{
			if (context.bbox[i][0][1] < 940.0 && context.bbox[i][0][1] > 130.0 && context.bbox[i][0][2] < 550.0)
			{
				continue;
			}
		}
		else 
		{
			if (context.bbox[i][0][1] < 840.0 && context.bbox[i][0][1] > 250.0 && context.bbox[i][0][2] > 0.0)
			{
				continue;
			}
		}

		if (context.bbox[i][tmp][2] > 350.0)
		{
			if (context.bbox[i][tmp][1] < 940.0 && context.bbox[i][tmp][1] > 130.0 && context.bbox[i][tmp][2] < 550.0)
			{
				continue;
			}
		}
		else
		{
			if (context.bbox[i][tmp][1] < 840.0 && context.bbox[i][tmp][1] > 250.0 && context.bbox[i][tmp][2] > 0.0)
			{
				continue;
			}
		}*/


		/*if ((context.bbox[i][0][1] < 560.0 && context.bbox[i][0][1] > 40.0 && context.bbox[i][0][2] < 260.0 && context.bbox[i][0][2] > 40.0) || (
			context.bbox[i][tmp][1] < 560.0 && context.bbox[i][tmp][1] > 40.0 && context.bbox[i][tmp][2] < 260.0 && context.bbox[i][tmp][2] > 40.0))
		{
			continue;
		}*/

		/*if (int(context.all_segs[i].size()) > 0)
		{
			visualize_state[i] = true;
			++visualize_tube_num;
		}*/ //------------------------------------------
		/*if ((context.bbox[i][0][1] < 770.0 && context.bbox[i][0][1] > 160.0 && context.bbox[i][0][2] < 530.0 && context.bbox[i][0][2] > 170.0) || (
			context.bbox[i][tmp][1] < 770.0 && context.bbox[i][tmp][1] > 160.0 && context.bbox[i][tmp][2] < 530.0 && context.bbox[i][tmp][2] > 170.0))
		{
			continue;
		}*/

		/*if ((context.bbox[i][0][1] < 800.0 && context.bbox[i][0][1] > 250.0 && context.bbox[i][0][2] < 530.0 && context.bbox[i][0][2] > 250.0) || (
			context.bbox[i][tmp][1] < 800.0 && context.bbox[i][tmp][1] > 250.0 && context.bbox[i][tmp][2] < 530.0 && context.bbox[i][tmp][2] > 250.0))
		{
			continue;
		}*/
		if (!context.visualize_state[i])
		{
			continue;
		}
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			//std::cout << i << "," << j << std::endl;

			int aa = context.all_segs[i][j]->aa_;
			int bb = context.all_segs[i][j]->bb_;
			double scale = context.all_segs[i][j]->scale_;


			for (int kk = aa; kk <= bb; kk++)
			{
				if (kk >= 0 && kk < synlen)
				{
					/* remain speed */

					//int k = context.all_segs[i][j]->interp_frame_No(kk);
					//cv::Rect roi = get_bbox_of_obj_by_frame(context.bbox, i, k);
					//VisMeta vm;
					//vm.src_frame_id = k; //原视频中的帧
					//vm.dst_frame_id = kk; //摘要视频中的帧
					//vm.scale = scale;
					//vm.box = roi;
					//vm.object_id = i + 1;
					//metas[k].push_back(vm);


					/* change speed */

					int k = context.all_segs[i][j]->interp_frame_No(kk);
					cv::Rect roi = get_bbox_of_obj_by_frame(context.bbox, i, k);
					int src_frame_idx = int(double((kk - aa)) / speed) + int(context.bbox[i][0][5]);
					VisMeta vm;
					vm.src_frame_id = src_frame_idx; //原视频中的帧
					vm.dst_frame_id = kk; //摘要视频中的帧
					vm.scale = scale;
					vm.box = roi;
					vm.object_id = i + 1;
					metas[src_frame_idx].push_back(vm);



					//cv::rectangle(frames[kk], roi, cv::Scalar(0, 255, 255));
					/*char text[1024];
					sprintf_s(text, 1024, "%d", i);
					cv::putText(frames[kk], text, cv::Point(roi.x + roi.width / 4, roi.y - 20), cv::FONT_HERSHEY_PLAIN, 1, cv::Scalar(0, 255, 255));*/
				}
			}
		}
	}
	//cap.set(cv::CAP_PROP_POS_FRAMES, 0);

	// 每帧中的bbox
	//std::map<int, std::vector<pair<int, cv::Vec6d>>>bbox_per_frame;
	// 前景像素数
	double fore_pix_num = 0.0;
	
	for (int epoch = 0; epoch < synlen; epoch += 1000)
	{
		int dst_start_frame_idx = epoch, dst_end_frame_idx = min(synlen - 1,epoch + 1000);

		std::vector<cv::Mat> frames;
		for (int i = 1; i <= dst_end_frame_idx-dst_start_frame_idx+1 /*context.synoplen*/; i++)
		{
			frames.push_back(tmp2.clone());
		}

		cv::VideoCapture cap_tmp(context.filepath + "/src.mp4");
		for (int i = 0; i < src_video_len; i++) //按照原视频中的帧
		{
			cv::Mat iimg;
			cap_tmp.read(iimg);

			if (i % 1000 == 0)
			{
				std::cout << i << "/" << src_video_len << std::endl;
			}

			if (metas[i].empty())
			{
				continue;
			}

			cv::resize(iimg, iimg, cv::Size(), context.scale, context.scale);

			//std::cout << iimg.rows << iimg.cols << std::endl;

			for (size_t j = 0; j < metas[i].size(); j++)
			{
				VisMeta vm = metas[i][j];
				cv::Rect roi = vm.box;
				int dx = int(roi.width * (vm.scale - 1) / 2.0);
				int dy = int(roi.height * (vm.scale - 1) / 2.0);

				cv::Rect new_roi(roi.x - dx, roi.y - dy, roi.width + 2 * dx, roi.height + 2 * dy);

				if (roi.width <= 0 || roi.height <= 0 || roi.width > video_width || roi.height > video_height)
					continue;

				if (new_roi.width <= 0 || new_roi.height <= 0 || new_roi.width > video_width || new_roi.height > video_height)
					continue;

				cv::Mat mask(new_roi.height, new_roi.width, CV_8U, cv::Scalar(1));
				cv::Mat new_img;
				cv::resize(iimg(roi), new_img, cv::Size(new_roi.width, new_roi.height));//iimg为原视频中的帧

				/*for (int u = 0; u < roi.height; u++)
				{
					for (int v = 0; v < roi.width; v++)
					{
						int uu = u + roi.y;
						int vv = v + roi.x;

						cv::Scalar csb = tmp2.at<cv::Vec3b>(uu, vv);
						cv::Scalar csi = iimg.at<cv::Vec3b>(uu, vv);
						double gb = csb[0] * 0.11 + csb[1] * 0.59 + csb[2] * 0.3;
						double gi = csi[0] * 0.11 + csi[1] * 0.59 + csi[2] * 0.3;

						if (fabs(gb - gi) >= 25)
						{
							mask.at<uchar>(u, v) = 1;
						}
					}
				}*/

				// 绘制patch
				if (vm.dst_frame_id >= dst_start_frame_idx && vm.dst_frame_id <= dst_end_frame_idx)
				{
					new_img.copyTo(frames[vm.dst_frame_id-dst_start_frame_idx](new_roi), mask); //dst摘要视频
				}
				

				//ID
				//cv::rectangle(frames[vm.dst_frame_id], roi, cv::Scalar(0, 0, 255));
				//char text[1024];
				//sprintf_s(text, 1024, "%d", i);
				//cv::putText(frames[vm.dst_frame_id], text, cv::Point(roi.x /*+ roi.width / 4*/, roi.y - 5), 
				//	cv::FONT_HERSHEY_PLAIN, 1, cv::Scalar(0, 0, 255));

			}

			for (size_t j = 0; j < metas[i].size(); j++)
			{
				VisMeta vm = metas[i][j];
				cv::Rect roi = vm.box;
				int dx = int(roi.width * (vm.scale - 1) / 2.0);
				int dy = int(roi.height * (vm.scale - 1) / 2.0);

				cv::Rect new_roi(roi.x - dx, roi.y - dy, roi.width + 2 * dx, roi.height + 2 * dy);

				if (roi.width <= 0 || roi.height <= 0 || roi.width > video_width || roi.height > video_height)
					continue;

				if (new_roi.width <= 0 || new_roi.height <= 0 || new_roi.width > video_width || new_roi.height > video_height)
					continue;

				//std::cout << roi << std::endl;

				//cv::Mat mask(new_roi.height, new_roi.width, CV_8U, cv::Scalar(1));

				/*cv::Mat new_img;
				cv::resize(iimg(roi), new_img, cv::Size(new_roi.width, new_roi.height));*/


				//new_img.copyTo(frames[vm.dst_frame_id](new_roi), mask);

				
				if (vm.dst_frame_id >= dst_start_frame_idx && vm.dst_frame_id <= dst_end_frame_idx)
				{
					// 绘制bbox
					//cv::Scalar color = random_color(tubes_group_id[vm.object_id] + 0 + 0 + 0); //------

					/*高亮occlusion*/
					cv::Scalar color = cv::Scalar(255, 255, 0);
					int line_thick = 1;
					// 高亮某一个tube
					/*if (vm.object_id == 11)
					{
						color = cv::Scalar(0, 255, 255);
						line_thick = 2;
						for (size_t seg_idx = 0; seg_idx < context.orig_colli[vm.object_id].size(); ++seg_idx)
						{
							if (context.orig_colli[vm.object_id][seg_idx].y <= vm.src_frame_id &&
								context.orig_colli[vm.object_id][seg_idx].z >= vm.src_frame_id)
							{
								color = cv::Scalar(0, 0, 255);
								break;
							}
						}
					}
					else {
						for (size_t seg_idx = 0; seg_idx < context.orig_colli[vm.object_id].size(); ++seg_idx)
						{
							if (context.orig_colli[vm.object_id][seg_idx].x == 11 && 
								context.orig_colli[vm.object_id][seg_idx].y <= vm.src_frame_id &&
								context.orig_colli[vm.object_id][seg_idx].z >= vm.src_frame_id)
							{
								color = cv::Scalar(0, 0, 255);
								break;
							}
						}
					}*/

					// 高亮occlusion
					for (size_t seg_idx = 0; seg_idx < context.orig_colli[vm.object_id].size(); ++seg_idx)
					{
						if (
							context.orig_colli[vm.object_id][seg_idx].y <= vm.src_frame_id &&
							context.orig_colli[vm.object_id][seg_idx].z >= vm.src_frame_id)
						{
							color = cv::Scalar(0, 0, 255);
							break;
						}
					}
					
					
					
					//cv::Scalar color = cv::Scalar(0, 0, 255);

					cv::rectangle(frames[vm.dst_frame_id-dst_start_frame_idx], new_roi, color, line_thick); //---------- bounding box
					char text[1024];

					// 打码
					cv::Scalar mask_color(160, 160, 160);
					int mask_thickness = -1;

					// src_sec
					//if (vm.object_id != 21 && vm.object_id != 44 && vm.object_id != 25 && vm.object_id != 41 && vm.object_id != 28)
					//{
					//	if (vm.object_id == 14 || vm.object_id == 24 || vm.object_id == 27 || vm.object_id == 37)
					//	{
					//		cv::Point center1(new_roi.x + new_roi.width / 4, new_roi.y + new_roi.height / 7);
					//		cv::Point center2(new_roi.x + 3 * new_roi.width / 4, new_roi.y + new_roi.height / 7);
					//		int radius = new_roi.width / 4;
					//		cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center1, radius, mask_color, mask_thickness);
					//		cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center2, radius, mask_color, mask_thickness);
					//	}
					//	else if (vm.object_id == 4)
					//	{
					//		if (vm.dst_frame_id - dst_start_frame_idx < 221) // our: 270, FSF:169, pccva: 1162, oct_speed: 221, oct_org: 258, input: 429
					//		{
					//			cv::Point center1(new_roi.x + 3 * new_roi.width / 4, new_roi.y + new_roi.height / 7);
					//			int radius = new_roi.width / 3.7;
					//			cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center1, radius, mask_color, mask_thickness);
					//		}
					//		else {
					//			cv::Point center1(new_roi.x + new_roi.width / 2, new_roi.y + new_roi.height / 7);
					//			int radius = new_roi.width / 3.7;
					//			cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center1, radius, mask_color, mask_thickness);
					//		}
					//		if (vm.dst_frame_id - dst_start_frame_idx < 211) // our: 260, FSF:159, pccva: 1152, oct_speed: 211, oct_org: 231, input: 402
					//		{
					//			cv::Rect rectangle2(new_roi.x, new_roi.y + new_roi.height / 5, new_roi.width / 2, new_roi.height / 4.8);
					//			cv::rectangle(frames[vm.dst_frame_id - dst_start_frame_idx], rectangle2, mask_color, mask_thickness);
					//		}
					//		else {
					//			cv::Rect rectangle2(new_roi.x + new_roi.width / 2, new_roi.y + new_roi.height / 3, new_roi.width / 2, new_roi.height / 5.5);
					//			cv::rectangle(frames[vm.dst_frame_id - dst_start_frame_idx], rectangle2, mask_color, mask_thickness);
					//		}
					//	}
					//	else {
					//		cv::Point center(new_roi.x + new_roi.width / 2, new_roi.y + new_roi.height / 7);
					//		int radius = new_roi.width / 3.5;
					//		cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center, radius, mask_color, mask_thickness);
					//	}
					//}

					// src_fir_1
					if (vm.object_id != 19 && vm.object_id != 28 && vm.object_id != 36)
					{
						if (vm.object_id != 20)
						{
							cv::Point center(new_roi.x + new_roi.width / 2, new_roi.y + new_roi.height / 7);
							int radius = new_roi.width / 4;
							cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center, radius, mask_color, mask_thickness);
						}
						else {
							cv::Point center(new_roi.x + new_roi.width / 4, new_roi.y + new_roi.height / 5);
							int radius = new_roi.width / 4;
							cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center, radius, mask_color, mask_thickness);
						}
					}

					//sprintf_s(text, 1024, "%d %d", vm.object_id + 0, tubes_group_id[vm.object_id] + 0 + 0 + 0); //------
					sprintf_s(text, 1024, "%d", vm.object_id + 0); //------

					cv::putText(frames[vm.dst_frame_id - dst_start_frame_idx], text, cv::Point(new_roi.x /*+ roi.width / 4*/, new_roi.y - 5),cv::FONT_HERSHEY_PLAIN, 1, color); //---------- tube number

					// 计算前景像素数
					fore_pix_num += new_roi.width * new_roi.height;

					cv::Vec6d bbox;
					bbox[0] = vm.dst_frame_id;
					bbox[1] = new_roi.x;
					bbox[2] = new_roi.y;
					bbox[3] = new_roi.width;
					bbox[4] = new_roi.height;
					bbox[5] = 0.0;

					context.bbox_per_frame[vm.dst_frame_id].push_back({vm.object_id, bbox});
				}
			}

		}

		for (size_t i = 0; i < frames.size(); i++)
		{
			if (i % 1000 == 0)
				std::cout << i << std::endl;
			capw.write(frames[i]);
		}
	}

	// 计算碰撞
	/*double colli_area = 0.0;
	for (auto it : bbox_per_frame)
	{
		for (int i = 0; i < int(it.second.size()) - 1; ++i)
		{
			for (int j = i + 1; j < int(it.second.size()); ++j)
			{
				if (tubes_group_id[it.second[i].first] != tubes_group_id[it.second[j].first])
				{
					colli_area += cal_colli_area(it.second[i].second, it.second[j].second);
				}
			}
		}
	}*/

	std::cout << "compute CA." << endl; // only for compared methods //--------------------------------------------------------
	std::map<int, std::map<int, int>> colli_relation;
	double colli_area = 0.0;
	for (auto it : context.bbox_per_frame)
	{
		for (int i = 0; i < int(it.second.size()) - 1; ++i)
		{
			for (int j = i + 1; j < int(it.second.size()); ++j)
			{
				if (tubes_group_id[it.second[i].first] != tubes_group_id[it.second[j].first])
				{
					double frame_colli = cal_colli_area(it.second[i].second, it.second[j].second);
					colli_area += frame_colli;
					if (frame_colli > 0)
					{
						if (colli_relation[it.second[i].first].count(it.second[j].first) == 0)
						{
							colli_relation[it.second[i].first][it.second[j].first] = 1;
						}
						else
						{
							colli_relation[it.second[i].first][it.second[j].first] += 1;
						}

						if (colli_relation[it.second[j].first].count(it.second[i].first) == 0)
						{
							colli_relation[it.second[j].first][it.second[i].first] = 1;
						}
						else
						{
							colli_relation[it.second[j].first][it.second[i].first] += 1;
						}
					}
				}
			}
		}
	}

	int colli_num = 0;
	for (auto it : colli_relation)
	{
		for (auto it2 : it.second)
		{
			if (it2.second > 30)
			{
				colli_num += 1;
			}
		}
	}
	std::cout << "colli_num: " << colli_num << std::endl;

	colli_num /= 2;

	std::cout << "colli_num: " << colli_num << std::endl;

	// 统计时序损失
	//int chrono_cost = 0;
	//for (int i = 0; i < int(context.all_segs.size()) - 1; i++)
	//{
	//	if (int(context.all_segs[i].size()) > 0 && context.visualize_state[i])
	//	{
	//		for (size_t j = i + 1; j < context.all_segs.size(); j++)
	//		{
	//			if (int(context.all_segs[j].size()) > 0 && context.visualize_state[j])
	//			{
	//				Segment* s1 = context.all_segs[i][0];
	//				Segment* s2 = context.all_segs[j][0];
	//				int d1 = s1->a_ - s2->a_;
	//				int d2 = s1->aa_ - s2->aa_;
	//				if (d1 * d2 < 0 && abs(d2) > chrono_slack)
	//				{
	//					++chrono_cost;
	//				}
	//				//ret += std::max(abs(abs(d1) - abs(d2)) - context.beg_slack, 0) * 100.0;
	//			}
	//		}
	//	}
	//}

	std::ofstream ofile_forePixNum(context.filepath + "/metrics.txt", ofstream::app);
	//ofile_forePixNum << "Collision: " << colli_area << "\n";
	ofile_forePixNum << "fore_pix_num: " << fore_pix_num << "\n";
	//ofile_forePixNum << "chrono_cost: " << chrono_cost << "\n";
	ofile_forePixNum << "visualize_tube_num: " << visualize_tube_num << "\n";
	ofile_forePixNum << "colli_num: " << colli_num << "\n";
	/*for (size_t i = 0; i < visualize_state.size(); ++i)
	{
		if (!visualize_state[i])
		{
			ofile_forePixNum << i + 1 << ",";
		}
		
	}*/

}


void visualize_for_ours(Context& context, int chrono_slack)
{
	int synlen = context.synoplen;

	cv::VideoCapture cap(context.filepath + "/src.mp4");

	int src_video_len = int(cap.get(cv::CAP_PROP_FRAME_COUNT));
	double fps = cap.get(cv::CAP_PROP_FPS);

	cv::Mat tmp2 = cv::imread(context.filepath + "/background.png");
	cv::resize(tmp2, tmp2, cv::Size(), context.scale, context.scale);
	int video_height = tmp2.rows;
	int video_width = tmp2.cols;

	cv::VideoWriter capw;

	if (!capw.open(context.filepath + "/synopsis.mp4", cv::VideoWriter::fourcc('m', 'p', 'v', '4'), fps, tmp2.size()))
	{
		std::cout << "Can not write!" << std::endl;
		return;
	}

	std::vector<std::vector<VisMeta>> metas;
	metas.resize(src_video_len);

	int visualize_tube_num = 0;
	//std::vector<bool>visualize_state(int(context.all_segs.size()), true);
	context.visualize_state.resize(int(context.all_segs.size()), true);



	
	//std::vector<int>vis_tube_vec = { 2, 4, 8, 10,11,14,15,17,18,19,20 };
	/*for (size_t i = 0; i < vis_tube_vec.size(); ++i)
	{
		context.visualize_state[vis_tube_vec[i] - 1] = true;
	}*/

	std::set<int>not_vis_tubes;

	

	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		// 不可视化边缘出现 边缘消失的tube和移动距离较小的tube
		/*double xmin = 1000000000000.0, xmax = -1.0, ymin = 100000000000000.0, ymax = -1.0;
		for (size_t ij = 0; ij < context.bbox[i].size(); ++ij)
		{
			xmin = min(xmin, context.bbox[i][ij][1]);
			xmax = max(xmax, context.bbox[i][ij][1]);

			ymin = min(ymin, context.bbox[i][ij][2]);
			ymax = max(ymax, context.bbox[i][ij][2]);
		}

		if ((xmax - xmin) < 140.0 && (ymax - ymin) < 140.0)
		{
			not_vis_tubes.insert(i);
			context.visualize_state[i] = false;
			continue;
		}*/

		int tmp = int(context.bbox[i].size()) - 1;
		// for fir1
		/*if (context.bbox[i][0][2] > 290.0)
		{
			if (context.bbox[i][0][1] < 770.0 && context.bbox[i][0][1] > 160.0 && context.bbox[i][0][2] < 530.0)
			{
				not_vis_tubes.insert(i);
				context.visualize_state[i] = false;

				
				continue;
			}
		}
		else
		{
			if (context.bbox[i][0][1] < 750.0 && context.bbox[i][0][1] > 160.0 && context.bbox[i][0][2] > 170.0)
			{
				not_vis_tubes.insert(i);
				context.visualize_state[i] = false;

				
				continue;
			}
		}

		if (context.bbox[i][tmp][2] > 290.0)
		{
			if (context.bbox[i][tmp][1] < 770.0 && context.bbox[i][tmp][1] > 160.0 && context.bbox[i][tmp][2] < 530.0)
			{
				not_vis_tubes.insert(i);
				context.visualize_state[i] = false;

				
				continue;
			}
		}
		else
		{
			if (context.bbox[i][tmp][1] < 750.0 && context.bbox[i][tmp][1] > 160.0 && context.bbox[i][tmp][2] > 170.0)
			{
				not_vis_tubes.insert(i);
				context.visualize_state[i] = false;

				
				continue;
			}
		}*/

		// for fir2
		/*if (context.bbox[i][0][2] > 350.0)
		{
			if (context.bbox[i][0][1] < 940.0 && context.bbox[i][0][1] > 130.0 && context.bbox[i][0][2] < 550.0)
			{
				not_vis_tubes.insert(i);
				context.visualize_state[i] = false;
				continue;
			}
		}
		else
		{
			if (context.bbox[i][0][1] < 840.0 && context.bbox[i][0][1] > 250.0 && context.bbox[i][0][2] > 0.0)
			{
				not_vis_tubes.insert(i);
				context.visualize_state[i] = false;
				continue;
			}
		}

		if (context.bbox[i][tmp][2] > 350.0)
		{
			if (context.bbox[i][tmp][1] < 940.0 && context.bbox[i][tmp][1] > 130.0 && context.bbox[i][tmp][2] < 550.0)
			{
				not_vis_tubes.insert(i);
				context.visualize_state[i] = false;
				continue;
			}
		}
		else
		{
			if (context.bbox[i][tmp][1] < 840.0 && context.bbox[i][tmp][1] > 250.0 && context.bbox[i][tmp][2] > 0.0)
			{
				not_vis_tubes.insert(i);
				context.visualize_state[i] = false;
				continue;
			}
		}*/

		// for sec
		/*if ((context.bbox[i][0][1] < 560.0 && context.bbox[i][0][1] > 40.0 && context.bbox[i][0][2] < 260.0 && context.bbox[i][0][2] > 40.0) || (
			context.bbox[i][tmp][1] < 560.0 && context.bbox[i][tmp][1] > 40.0 && context.bbox[i][tmp][2] < 260.0 && context.bbox[i][tmp][2] > 40.0))
		{
			continue;
		}*/

		//

		// old
		/*if ((context.bbox[i][0][1] < 770.0 && context.bbox[i][0][1] > 160.0 && context.bbox[i][0][2] < 530.0 && context.bbox[i][0][2] > 170.0) || (
			context.bbox[i][tmp][1] < 770.0 && context.bbox[i][tmp][1] > 160.0 && context.bbox[i][tmp][2] < 530.0 && context.bbox[i][tmp][2] > 170.0))
		{
			continue;
		}*/

		/*if ((context.bbox[i][0][1] < 800.0 && context.bbox[i][0][1] > 250.0 && context.bbox[i][0][2] < 530.0 && context.bbox[i][0][2] > 250.0) || (
			context.bbox[i][tmp][1] < 800.0 && context.bbox[i][tmp][1] > 250.0 && context.bbox[i][tmp][2] < 530.0 && context.bbox[i][tmp][2] > 250.0))
		{
			continue;
		}*/
	}

	// 不可视化与已经确定不可视化的tube有occlusion的tube
	/*int tube_beg = 115;
	std::vector<int>not_vis_num = { 122,123, 134 ,136 ,157 ,180 ,184 ,203 ,204,140,207,130,135,177,
		191,179,143,192,202,206,147,209,208,
	200};*/

	/*for (size_t i = 0; i < not_vis_num.size(); ++i)
	{
		context.visualize_state[not_vis_num[i]-tube_beg] = false;
		not_vis_tubes.insert(not_vis_num[i]-tube_beg);
	}*/
	
	
	/*while (!not_vis_tubes.empty())
	{
		std::set<int>not_vis_tubes_next;
		for (std::set<int>::iterator iter = not_vis_tubes.begin(); iter != not_vis_tubes.end(); ++iter)
		{
			int tube_id = *iter;
			for (size_t i = 0; i < context.all_segs[tube_id].size(); ++i)
			{
				for (size_t j = 0; j < context.all_segs[tube_id][i]->colli_segs_.size(); ++j)
				{
					int colli_tube_id = context.all_segs[tube_id][i]->colli_segs_[j]->tube_id_;
					if (context.visualize_state[colli_tube_id])
					{
						context.visualize_state[colli_tube_id] = false;
						not_vis_tubes_next.insert(colli_tube_id);
					}
				}
			}
		}
		not_vis_tubes = not_vis_tubes_next;
	}*/

	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		/*if (context.circulated_tube[i] == 1)
		{
			continue;
		}*/
		if (!context.visualize_state[i])
		{
			continue;
		}
		++visualize_tube_num;

		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			//std::cout << i << "," << j << std::endl;

			int aa = context.all_segs[i][j]->aa_;
			int bb = context.all_segs[i][j]->bb_;
			double scale = context.all_segs[i][j]->scale_;

			for (int kk = aa; kk <= bb; kk++)
			{
				if (kk >= 0 && kk < synlen)
				{
					int k = context.all_segs[i][j]->interp_frame_No(kk);
					cv::Rect roi = get_bbox_of_obj_by_frame(context.bbox, i, k);
					
					//std::cout << src_frame_idx << std::endl;
					/*cap.set(cv::CAP_PROP_POS_FRAMES, k);
					cv::Mat iimg;
					cap >> iimg;
					cv::resize(iimg, iimg, cv::Size(), context.scale, context.scale);

					cv::Mat mask(roi.height, roi.width, CV_8U, cv::Scalar(1));*/

					/*for (int u = 0; u < roi.height; u++)
					{
						for (int v = 0; v < roi.width; v++)
						{
							int uu = u + roi.y;
							int vv = v + roi.x;

							cv::Scalar csb = tmp2.at<cv::Vec3b>(uu, vv);
							cv::Scalar csi = iimg.at<cv::Vec3b>(uu, vv);
							double gb = csb[0] * 0.11 + csb[1] * 0.59 + csb[2] * 0.3;
							double gi = csi[0] * 0.11 + csi[1] * 0.59 + csi[2] * 0.3;

							if (fabs(gb - gi) >= 25)
							{
								mask.at<uchar>(u, v) = 1;
							}
						}
					}*/

					//iimg(roi).copyTo(frames[kk](roi), mask);

					VisMeta vm;
					vm.src_frame_id = k; //原视频中的帧
					vm.dst_frame_id = kk; //摘要视频中的帧
					vm.scale = scale;
					vm.box = roi;
					vm.object_id = i + 1;
					metas[k].push_back(vm);

					//cv::rectangle(frames[kk], roi, cv::Scalar(0, 255, 255));
					/*char text[1024];
					sprintf_s(text, 1024, "%d", i);
					cv::putText(frames[kk], text, cv::Point(roi.x + roi.width / 4, roi.y - 20), cv::FONT_HERSHEY_PLAIN, 1, cv::Scalar(0, 255, 255));*/
				}
			}
		}
	}
	//cap.set(cv::CAP_PROP_POS_FRAMES, 0);

	// 每帧中的bbox
	//std::map<int, std::vector<pair<int, cv::Vec6d>>>bbox_per_frame;
	// 前景像素数
	double fore_pix_num = 0.0;
	for (int epoch = 0; epoch < synlen; epoch += 1000)
	{
		int dst_start_frame_idx = epoch, dst_end_frame_idx = min(synlen - 1, epoch + 1000);

		std::vector<cv::Mat> frames;
		for (int i = 1; i <= dst_end_frame_idx - dst_start_frame_idx + 1; i++)
		{
			frames.push_back(tmp2.clone());
		}

		cv::VideoCapture cap_tmp(context.filepath + "/src.mp4");
		for (int i = 0; i < src_video_len; i++) //按照原视频中的帧
		{
			cv::Mat iimg;
			cap_tmp.read(iimg);

			if (i % 1000 == 0)
			{
				std::cout << i << "/" << src_video_len << std::endl;
			}

			if (metas[i].empty())
			{
				continue;
			}

			cv::resize(iimg, iimg, cv::Size(), context.scale, context.scale);

			//std::cout << iimg.rows << iimg.cols << std::endl;

			for (size_t j = 0; j < metas[i].size(); j++)
			{
				VisMeta vm = metas[i][j];
				cv::Rect roi = vm.box;
				int dx = int(roi.width * (vm.scale - 1) / 2.0);
				int dy = int(roi.height * (vm.scale - 1) / 2.0);

				cv::Rect new_roi(roi.x - dx, roi.y - dy, roi.width + 2 * dx, roi.height + 2 * dy);

				if (roi.width <= 0 || roi.height <= 0 || roi.width > video_width || roi.height > video_height)
					continue;

				if (new_roi.width <= 0 || new_roi.height <= 0 || new_roi.width > video_width || new_roi.height > video_height)
					continue;

				cv::Mat mask(new_roi.height, new_roi.width, CV_8U, cv::Scalar(1));
				cv::Mat new_img;
				cv::resize(iimg(roi), new_img, cv::Size(new_roi.width, new_roi.height));//iimg为原视频中的帧

				/*for (int u = 0; u < roi.height; u++)
				{
					for (int v = 0; v < roi.width; v++)
					{
						int uu = u + roi.y;
						int vv = v + roi.x;

						cv::Scalar csb = tmp2.at<cv::Vec3b>(uu, vv);
						cv::Scalar csi = iimg.at<cv::Vec3b>(uu, vv);
						double gb = csb[0] * 0.11 + csb[1] * 0.59 + csb[2] * 0.3;
						double gi = csi[0] * 0.11 + csi[1] * 0.59 + csi[2] * 0.3;

						if (fabs(gb - gi) >= 25)
						{
							mask.at<uchar>(u, v) = 1;
						}
					}
				}*/

				// 绘制patch
				if (vm.dst_frame_id >= dst_start_frame_idx && vm.dst_frame_id <= dst_end_frame_idx)
				{
					new_img.copyTo(frames[vm.dst_frame_id - dst_start_frame_idx](new_roi), mask); //dst摘要视频
				}


				//ID
				//cv::rectangle(frames[vm.dst_frame_id], roi, cv::Scalar(0, 0, 255));
				//char text[1024];
				//sprintf_s(text, 1024, "%d", i);
				//cv::putText(frames[vm.dst_frame_id], text, cv::Point(roi.x /*+ roi.width / 4*/, roi.y - 5), 
				//	cv::FONT_HERSHEY_PLAIN, 1, cv::Scalar(0, 0, 255));

			}

			for (size_t j = 0; j < metas[i].size(); j++)
			{
				VisMeta vm = metas[i][j];
				cv::Rect roi = vm.box;
				int dx = int(roi.width * (vm.scale - 1) / 2.0);
				int dy = int(roi.height * (vm.scale - 1) / 2.0);

				cv::Rect new_roi(roi.x - dx, roi.y - dy, roi.width + 2 * dx, roi.height + 2 * dy);

				if (roi.width <= 0 || roi.height <= 0 || roi.width > video_width || roi.height > video_height)
					continue;

				if (new_roi.width <= 0 || new_roi.height <= 0 || new_roi.width > video_width || new_roi.height > video_height)
					continue;

				//std::cout << roi << std::endl;

				//cv::Mat mask(new_roi.height, new_roi.width, CV_8U, cv::Scalar(1));

				/*cv::Mat new_img;
				cv::resize(iimg(roi), new_img, cv::Size(new_roi.width, new_roi.height));*/


				//new_img.copyTo(frames[vm.dst_frame_id](new_roi), mask);


				if (vm.dst_frame_id >= dst_start_frame_idx && vm.dst_frame_id <= dst_end_frame_idx)
				{
					// 绘制bbox
					//cv::Scalar color = random_color(tubes_group_id[vm.object_id] + 0 + 0 + 0); //------

					/*高亮occlusion*/
					cv::Scalar color = cv::Scalar(255, 255, 0);
					int line_thick = 1;
					// 高亮某一个tube
					/*if (vm.object_id == 11)
					{
						color = cv::Scalar(0, 255, 255);
						line_thick = 2;
						for (size_t seg_idx = 0; seg_idx < context.orig_colli[vm.object_id].size(); ++seg_idx)
						{
							if (context.orig_colli[vm.object_id][seg_idx].y <= vm.src_frame_id &&
								context.orig_colli[vm.object_id][seg_idx].z >= vm.src_frame_id)
							{
								color = cv::Scalar(0, 0, 255);
								break;
							}
						}
					}
					else {
						for (size_t seg_idx = 0; seg_idx < context.orig_colli[vm.object_id].size(); ++seg_idx)
						{
							if (context.orig_colli[vm.object_id][seg_idx].x == 11 &&
								context.orig_colli[vm.object_id][seg_idx].y <= vm.src_frame_id &&
								context.orig_colli[vm.object_id][seg_idx].z >= vm.src_frame_id)
							{
								color = cv::Scalar(0, 0, 255);
								break;
							}
						}
					}*/

					// 高亮occlusion
					for (size_t seg_idx = 0; seg_idx < context.orig_colli[vm.object_id].size(); ++seg_idx)
					{
						if (
							context.orig_colli[vm.object_id][seg_idx].y <= vm.src_frame_id &&
							context.orig_colli[vm.object_id][seg_idx].z >= vm.src_frame_id)
						{
							color = cv::Scalar(0, 0, 255);
							break;
						}
					}



					//cv::Scalar color = cv::Scalar(0, 0, 255);

					cv::rectangle(frames[vm.dst_frame_id - dst_start_frame_idx], new_roi, color, line_thick);
					char text[1024];

					// 打码
					
					cv::Scalar mask_color(160, 160, 160);
					int mask_thickness = -1;
					
					// src_sec
					/*if (vm.object_id != 21 && vm.object_id != 44 && vm.object_id != 25 && vm.object_id != 41 && vm.object_id != 28)
					{
						if (vm.object_id == 14 || vm.object_id == 24 || vm.object_id == 27 || vm.object_id == 37)
						{
							cv::Point center1(new_roi.x + new_roi.width / 4, new_roi.y + new_roi.height / 7);
							cv::Point center2(new_roi.x + 3 * new_roi.width / 4, new_roi.y + new_roi.height / 7);
							int radius = new_roi.width / 4;
							cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center1, radius, mask_color, mask_thickness);
							cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center2, radius, mask_color, mask_thickness);
						}
						else if (vm.object_id == 4)
						{
							if (vm.dst_frame_id - dst_start_frame_idx < 270)
							{
								cv::Point center1(new_roi.x + 3 * new_roi.width / 4, new_roi.y + new_roi.height / 7);
								int radius = new_roi.width / 3.7;
								cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center1, radius, mask_color, mask_thickness);
							}
							else {
								cv::Point center1(new_roi.x + new_roi.width / 2, new_roi.y + new_roi.height / 7);
								
								int radius = new_roi.width / 3.7;
								cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center1, radius, mask_color, mask_thickness);	
							}
							if (vm.dst_frame_id - dst_start_frame_idx < 260)
							{
								cv::Rect rectangle2(new_roi.x, new_roi.y + new_roi.height / 5, new_roi.width / 2, new_roi.height / 4.8);
								cv::rectangle(frames[vm.dst_frame_id - dst_start_frame_idx], rectangle2, mask_color, mask_thickness);
							}
							else {
								cv::Rect rectangle2(new_roi.x + new_roi.width / 2, new_roi.y + new_roi.height / 3, new_roi.width / 2, new_roi.height / 5.5);
								cv::rectangle(frames[vm.dst_frame_id - dst_start_frame_idx], rectangle2, mask_color, mask_thickness);
							}
						}
						else {
							cv::Point center(new_roi.x + new_roi.width / 2, new_roi.y + new_roi.height / 7);
							int radius = new_roi.width / 3.5;
							cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center, radius, mask_color, mask_thickness);
						}
					}*/

					// src_fir_1
					if (vm.object_id != 19 && vm.object_id != 28 && vm.object_id != 36 && vm.object_id != 29)
					{
						if (vm.object_id != 20)
						{
							cv::Point center(new_roi.x + new_roi.width / 2, new_roi.y + new_roi.height / 7);
							int radius = new_roi.width / 4;
							cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center, radius, mask_color, mask_thickness);
						}
						else {
							cv::Point center(new_roi.x + new_roi.width / 4, new_roi.y + new_roi.height / 5);
							int radius = new_roi.width / 4;
							cv::circle(frames[vm.dst_frame_id - dst_start_frame_idx], center, radius, mask_color, mask_thickness);
						}
					}
					


					//sprintf_s(text, 1024, "%d %d", vm.object_id + 0, tubes_group_id[vm.object_id] + 0 + 0 + 0); //------
					sprintf_s(text, 1024, "%d", vm.object_id + 0); //------

					cv::putText(frames[vm.dst_frame_id - dst_start_frame_idx], text, cv::Point(new_roi.x /*+ roi.width / 4*/, new_roi.y - 5),
						cv::FONT_HERSHEY_PLAIN, 1, color);

					// 计算前景像素数
					fore_pix_num += new_roi.width * new_roi.height;

					cv::Vec6d bbox;
					bbox[0] = vm.dst_frame_id;
					bbox[1] = new_roi.x;
					bbox[2] = new_roi.y;
					bbox[3] = new_roi.width;
					bbox[4] = new_roi.height;
					bbox[5] = 0.0;

					context.bbox_per_frame[vm.dst_frame_id].push_back({ vm.object_id, bbox });
				}
			}

		}

		for (size_t i = 0; i < frames.size(); i++)
		{
			if (i % 1000 == 0)
				std::cout << i << std::endl;
			capw.write(frames[i]);
		}
	}

	// 计算视觉碰撞
	//double colli_area = compute_energy_collision4_vis(context);
	/*for (auto it : bbox_per_frame)
	{
		for (int i = 0; i < int(it.second.size()) - 1; ++i)
		{
			for (int j = i + 1; j < int(it.second.size()); ++j)
			{
				
				colli_area += cal_colli_area(it.second[i].second, it.second[j].second);
				
			}
		}
	}*/

	std::cout << "compute CA." << endl; // only for compared methods //--------------------------------------------------------
	std::map<int, std::map<int,int>> colli_relation;
	double colli_area = 0.0;
	for (auto it : context.bbox_per_frame)
	{
		for (int i = 0; i < int(it.second.size()) - 1; ++i)
		{
			for (int j = i + 1; j < int(it.second.size()); ++j)
			{
				if (context.occlu_graph[it.second[i].first].count(it.second[j].first) == 0)
				{
					double frame_colli = cal_colli_area(it.second[i].second, it.second[j].second);
					colli_area += frame_colli;
					if (frame_colli > 0)
					{
						if (colli_relation[it.second[i].first].count(it.second[j].first) == 0)
						{
							colli_relation[it.second[i].first][it.second[j].first] = 1;
						}
						else
						{
							colli_relation[it.second[i].first][it.second[j].first] += 1;
						}
						
						if (colli_relation[it.second[j].first].count(it.second[i].first) == 0)
						{
							colli_relation[it.second[j].first][it.second[i].first] = 1;
						}
						else
						{
							colli_relation[it.second[j].first][it.second[i].first] += 1;
						}
					}
				}
			}
		}
	}

	int colli_num = 0;
	for (auto it : colli_relation)
	{
		for (auto it2 : it.second)
		{
			if (it2.second > 30)
			{
				colli_num += 1;
			}
		}
	}
	std::cout << "colli_num: " << colli_num << std::endl;

	colli_num /= 2;

	std::cout << "colli_num: " << colli_num << std::endl;

	
	std::ofstream ofile_forePixNum(context.filepath + "/metrics.txt", ofstream::app);
	
	ofile_forePixNum << "fore_pix_num: "<< fore_pix_num << "\n";
	
	ofile_forePixNum << "visualize_tube_num: " << visualize_tube_num << "\n";

	ofile_forePixNum << "colli_num: " << colli_num << "\n";
}

void visualize_for_ours_step(Context& context)
{
	context.visualize_state.resize(int(context.all_segs.size()), true);
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		// 不可视化边缘出现 边缘消失的tube和移动距离较小的tube
		/*double xmin = 1000000000000.0, xmax = -1.0, ymin = 100000000000000.0, ymax = -1.0;
		for (size_t ij = 0; ij < context.bbox[i].size(); ++ij)
		{
			xmin = min(xmin, context.bbox[i][ij][1]);
			xmax = max(xmax, context.bbox[i][ij][1]);

			ymin = min(ymin, context.bbox[i][ij][2]);
			ymax = max(ymax, context.bbox[i][ij][2]);
		}

		if ((xmax - xmin) < 140.0 && (ymax - ymin) < 140.0)
		{
			context.visualize_state[i] = false;
			continue;
		}*/

		int tmp = int(context.bbox[i].size()) - 1;
		// for fir1
		/*if (context.bbox[i][0][2] > 290.0)
		{
			if (context.bbox[i][0][1] < 770.0 && context.bbox[i][0][1] > 160.0 && context.bbox[i][0][2] < 530.0)
			{
				context.visualize_state[i] = false;
				continue;
			}
		}
		else
		{
			if (context.bbox[i][0][1] < 750.0 && context.bbox[i][0][1] > 160.0 && context.bbox[i][0][2] > 170.0)
			{
				context.visualize_state[i] = false;
				continue;
			}
		}

		if (context.bbox[i][tmp][2] > 290.0)
		{
			if (context.bbox[i][tmp][1] < 770.0 && context.bbox[i][tmp][1] > 160.0 && context.bbox[i][tmp][2] < 530.0)
			{
				context.visualize_state[i] = false;
				continue;
			}
		}
		else
		{
			if (context.bbox[i][tmp][1] < 750.0 && context.bbox[i][tmp][1] > 160.0 && context.bbox[i][tmp][2] > 170.0)
			{
				context.visualize_state[i] = false;
				continue;
			}
		}*/

		// for fir2
		/*if (context.bbox[i][0][2] > 350.0)
		{
			if (context.bbox[i][0][1] < 940.0 && context.bbox[i][0][1] > 130.0 && context.bbox[i][0][2] < 550.0)
			{
				context.visualize_state[i] = false;
				continue;
			}
		}
		else
		{
			if (context.bbox[i][0][1] < 840.0 && context.bbox[i][0][1] > 250.0 && context.bbox[i][0][2] > 0.0)
			{
				context.visualize_state[i] = false;
				continue;
			}
		}

		if (context.bbox[i][tmp][2] > 350.0)
		{
			if (context.bbox[i][tmp][1] < 940.0 && context.bbox[i][tmp][1] > 130.0 && context.bbox[i][tmp][2] < 550.0)
			{
				context.visualize_state[i] = false;
				continue;
			}
		}
		else
		{
			if (context.bbox[i][tmp][1] < 840.0 && context.bbox[i][tmp][1] > 250.0 && context.bbox[i][tmp][2] > 0.0)
			{
				context.visualize_state[i] = false;
				continue;
			}
		}*/

		// for sec
		/*if ((context.bbox[i][0][1] < 560.0 && context.bbox[i][0][1] > 40.0 && context.bbox[i][0][2] < 260.0 && context.bbox[i][0][2] > 40.0) || (
			context.bbox[i][tmp][1] < 560.0 && context.bbox[i][tmp][1] > 40.0 && context.bbox[i][tmp][2] < 260.0 && context.bbox[i][tmp][2] > 40.0))
		{
			continue;
		}*/
	}


	// 统计时序损失
	int chrono_cost = 0;
	for (int i = 0; i < int(context.all_segs.size()) - 1; i++)
	{
		if (int(context.all_segs[i].size()) > 0 && context.visualize_state[i])
		{
			for (size_t j = i + 1; j < context.all_segs.size(); j++)
			{
				if (int(context.all_segs[j].size()) > 0 && context.visualize_state[j])
				{
					Segment* s1 = context.all_segs[i][0];
					Segment* s2 = context.all_segs[j][0];
					int d1 = s1->a_ - s2->a_;
					int d2 = s1->aa_ - s2->aa_;
					if (d1 * d2 < 0 && abs(d2) > 200)
					{
						++chrono_cost;
					}
				}
			}
		}
	}

	// 计算视觉碰撞
	double colli_area = compute_energy_collision4_vis(context);
	

	std::ofstream ofile_forePixNum(context.filepath + "/metrics.txt", ofstream::app);
	ofile_forePixNum << "mcmc2_iteration_num: " << context.mcmc2_iteration_num << "\n";
	ofile_forePixNum << "chrono_cost: " << chrono_cost << "\n";
	ofile_forePixNum << "vis collision: " << colli_area << "\n";
}

void visualize_root(string filepath, std::vector<Segment*>& segs, std::vector<std::vector<cv::Vec6d>>& bbox, double scale,
	int root_id, int video_len)
{
	cv::VideoCapture cap(filepath + "/src.mp4");

	int max_bb = -10000;
	for (size_t i = 0; i < segs.size(); i++)
	{
		if (segs[i]->bb_ > max_bb)
		{
			max_bb = segs[i]->bb_;
		}
	}

	std::cout << max_bb << std::endl;

	std::vector<cv::Mat> frames;
	cv::Mat tmp2 = cv::imread(filepath + "/background.png");
	cv::resize(tmp2, tmp2, cv::Size(), scale, scale);
	for (int i = 0; i < video_len; i++)
	{
		frames.push_back(tmp2.clone());
	}

	for (size_t i = 0; i < segs.size(); i++)
	{
		cout << "Seg:" << i << std::endl;
		int aa = segs[i]->aa_;
		int bb = segs[i]->bb_;

		cout << aa << "," << bb <<"++++++++++"<< std::endl;

		for (int kk = aa; kk <= bb; kk++)
		{
			int k = segs[i]->interp_frame_No(kk);

			std::cout << "   " << kk << "," << k << std::endl;
			cv::Rect roi = get_bbox_of_obj_by_frame(bbox, segs[i]->tube_id_, k);
			cap.set(cv::CAP_PROP_POS_FRAMES, k);
			cv::Mat iimg;
			cap >> iimg;
			cv::resize(iimg, iimg, cv::Size(), scale, scale);
			if (kk >= 0 && kk < video_len)
			{
				iimg(roi).copyTo(frames[kk](roi));
				//cv::rectangle(frames[kk], roi, cv::Scalar(255, 0, 0));
			}

		}
	}

	cv::VideoWriter capw;
	char filename[1024];
	sprintf_s(filename, 1024, "%s/tube_videos/%d.mp4", filepath.c_str(), root_id);
	std::cout << capw.open(filename, cv::VideoWriter::fourcc('H', '2', '6', '4'), 30, frames[0].size()) << std::endl;
	for (size_t i = 0; i < frames.size(); i++)
	{
		std::cout << i << std::endl;
		capw.write(frames[i]);
	}
}


void visualize_A_frame(cv::Mat & frame, int frame_id, Context & context)
{
	cv::Mat bg = cv::imread(context.filepath + "/background.png");
	frame = bg.clone();
	cv::VideoCapture cap(context.filepath + "/src.mp4");

	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{

			int aa = context.all_segs[i][j]->aa_;
			int bb = context.all_segs[i][j]->bb_;

			if (frame_id >= aa && frame_id <= bb)
			{
				int k = context.all_segs[i][j]->interp_frame_No(frame_id);
				cv::Rect roi = get_bbox_of_obj_by_frame(context.bbox, i, k);
				cap.set(cv::CAP_PROP_POS_FRAMES, k);
				cv::Mat iimg;
				cap >> iimg;
				

				cv::Mat mask(roi.height, roi.width, CV_8U, cv::Scalar(0));

				for (int u = 0; u < roi.height; u++)
				{
					for (int v = 0; v < roi.width; v++)
					{
						int uu = u + roi.y;
						int vv = v + roi.x;

						cv::Scalar csb = bg.at<cv::Vec3b>(uu, vv);
						cv::Scalar csi = iimg.at<cv::Vec3b>(uu, vv);
						double gb = csb[0] * 0.11 + csb[1] * 0.59 + csb[2] * 0.3;
						double gi = csi[0] * 0.11 + csi[1] * 0.59 + csi[2] * 0.3;

						if (fabs(gb - gi) >= 25)
						{
							mask.at<uchar>(u, v) = 1;
						}
					}
				}
				iimg(roi).copyTo(frame(roi), mask);

				//cv::rectangle(frames[kk], roi, cv::Scalar(0, 255, 255));
				char text[1024];
				sprintf_s(text, 1024, "%d", i);
				cv::putText(frame, text, cv::Point(roi.x + roi.width / 4, roi.y + 20), cv::FONT_HERSHEY_PLAIN, 1, cv::Scalar(255, 0, 0));
				
			}
		}
	}
}

