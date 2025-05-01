#include "fast_forward.h"
#include <opencv2/opencv.hpp>


void fast_forward(const std::string& src_video_path, const std::string& dst_video_path, int synlen, double scale)
{
	cv::VideoCapture cap(src_video_path);
	int src_video_len = int(cap.get(cv::CAP_PROP_FRAME_COUNT));
	double fps = cap.get(cv::CAP_PROP_FPS);
	int video_width = int(cap.get(cv::CAP_PROP_FRAME_WIDTH));
	int video_height = int(cap.get(cv::CAP_PROP_FRAME_HEIGHT));	
	int a = src_video_len / synlen;

	cv::VideoWriter capw;

	if (!capw.open(dst_video_path, cv::VideoWriter::fourcc('H', '2', '6', '4'), fps, 
		cv::Size(int(video_width * scale),int(video_height * scale))))
	{
		std::cout << "Can not write!" << std::endl;
		return;
	}

	int dst_video_width = int(capw.get(cv::CAP_PROP_FRAME_WIDTH));
	int dst_video_height = int(capw.get(cv::CAP_PROP_FRAME_HEIGHT));

	int b = 0;
	while (b < src_video_len)
	{
		cv::Mat iimg;
		cap.read(iimg);
		if (b % 1000 == 0)
		{
			std::cout << b << "/" << src_video_len << std::endl;
		}

		if (b % a == 0)
		{
			cv::resize(iimg, iimg, cv::Size(), scale, scale);

			/*std::cout << dst_video_width << "," << dst_video_width << std::endl;
			std::cout << iimg.cols << "," << iimg.rows << std::endl;*/

			capw.write(iimg);
		}
		b++;
	}
}