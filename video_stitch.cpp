#include "video_stitch.h"
#include <opencv2/opencv.hpp>

void stich_two_videos(const std::string& path1, const std::string& path2, const std::string& out_path)
{
	cv::VideoCapture cap1(path1);
	cv::VideoCapture cap2(path2);

	int src_video_len1 = int(cap1.get(cv::CAP_PROP_FRAME_COUNT));
	int src_video_len2 = int(cap2.get(cv::CAP_PROP_FRAME_COUNT));

	int dst_video_len = src_video_len1 < src_video_len2 ? src_video_len1: src_video_len2;

	double fps = cap1.get(cv::CAP_PROP_FPS);

	int w1 = int(cap1.get(cv::CAP_PROP_FRAME_WIDTH));
	int h1 = int(cap1.get(cv::CAP_PROP_FRAME_HEIGHT));
	int w2 = int(cap1.get(cv::CAP_PROP_FRAME_WIDTH));
	int h2 = int(cap1.get(cv::CAP_PROP_FRAME_HEIGHT));

	
	cv::Rect box1(0, 0, w1, h1);
	cv::Rect box2_from(0, 0, w2, h2);
	cv::Rect box2_to(w1 + 50, 0, w2, h2);

	int w = w1 + w2 + 100;
	int h = h1 + 100;
	
	cv::VideoWriter capw;
	capw.open(out_path, cv::VideoWriter::fourcc('H', '2', '6', '4'), fps,
		cv::Size(w, h));

	int a = 0;
	while (a < dst_video_len)
	{
		if (a % 1000 == 0)
		{
			std::cout << a << "/" << dst_video_len << std::endl;
		}

		cv::Mat img1, img2;
		cap1.read(img1);
		cap2.read(img2);

		cv::Mat img(h, w, CV_8UC3, cv::Scalar(0, 0, 0));

		img1(box1).copyTo(img(box1));
		img2(box2_from).copyTo(img(box2_to));

		capw.write(img);
		a++;
	}
}