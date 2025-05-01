#include "loadata.h"
#include <fstream>

void read_activity(string filepath, std::vector<std::vector<double>>& activity, int tube_num)
{
	char filename[1024];
	for (int i = 1; i <= tube_num; i++) // 读取所有value文件，每个value文件对应一个tube所有出现帧中的activity值
	{
		sprintf_s(filename, 1024, "%s/%dvalue.txt", filepath.c_str(), i);
		std::ifstream infile(filename);
		std::vector<double> tmp;
		std::string line;
		while (getline(infile, line) && line != "")
		{
			double a, b;
			std::istringstream(line) >> a >> b;
			tmp.push_back(b);
		}
		activity.push_back(tmp);
	}
}

void read_boundingbox(string filepath, std::vector<std::vector<cv::Vec6d>>& bbox, int tube_num, double scale, int video_height,
	int video_width)
{
	char filename[1024];
	for (int i = 1; i <= tube_num; i++) // 读取所有node文件，每个node文件对应一个tube所有出现帧中的bbox
	{
		sprintf_s(filename, 1024, "%s/%dnode.txt", filepath.c_str(), i);
		std::ifstream infile(filename);
		std::vector<cv::Vec6d> tmp;
		std::string line;
		while (getline(infile, line) && line != "")
		{
			double a, b, c, d, e;
			std::istringstream(line) >> a >> b >> c >> d >> e;
			
			/*if (b < 0)
				b = 0;
			if (b >= video_width)
				b = (double)video_width - 1;
			if (c < 0)
				c = 0;
			if (c >= video_height)
				c = (double)video_height - 1;
			
			double bb = b + d;
			double cc = c + e;
			if (bb < 0)
				bb = 0;
			if (cc < 0)
				cc = 0;
			if (bb >= video_width)
				bb = video_width - 1;
			if (cc >= video_height)
				cc = video_height - 1;

			d = bb - b + 1;
			e = cc - c + 1;

			if (d < 0)
				d = 0;
			if (e < 0)
				e = 0;*/


			cv::Vec6d bbb;
			bbb[0] = a;
			bbb[1] = b*scale;
			bbb[2] = c*scale;
			bbb[3] = d*scale;
			bbb[4] = e*scale;
			bbb[5] = 0;
			tmp.push_back(bbb);
		}
		bbox.push_back(tmp);
	}
}

void read_bestCopy(string filepath, std::vector<std::vector<Segment*>>& all_segs)
{
	string filename = filepath + "/bestCopy.txt";
	std::ifstream infile(filename);

	std::string line;
	int tube_id = -1;
	std::vector<Segment*>segs;
	int beg = 0, end = 0;
	while (getline(infile, line) && line != "")
	{
		// 一个新的tube
		if (int(line.length()) < 7)
		{
			if (tube_id >= 0)
			{
				std::vector<Segment*>tube_segs;
				tube_segs.assign(segs.begin() + beg, segs.begin() + end);
				all_segs.push_back(tube_segs);

				beg = end;
			}
			
			++tube_id;
		}
		else 
		{
			// 读这个tube的seg
			int sum = 0, i = 0;
			int num = -1, a = -1, b = -1, aa = -1, bb = -1;
			for (size_t i = 0; i < line.size(); ++i)
			{
				if (isdigit(line[i]))
				{
					sum = sum * 10 + line[i] - '0'; // 转整形
				}
				else
				{
					if (num == -1)
					{
						num = sum;
					}
					else if (a == -1)
					{
						a = sum;
					}
					else if (b == -1)
					{
						b = sum;
					}
					else if (aa == -1)
					{
						aa = sum;
					}
					else
					{
						bb = sum;
					}
					sum = 0;
				}
			}
			bb = sum;
			sum = 0;

			segs.push_back(new Segment(a, b, aa, bb, tube_id, bb - aa + 1));
			++end;
		}
	}

	std::vector<Segment*>tube_segs;
	tube_segs.assign(segs.begin() + beg, segs.begin() + end);
	all_segs.push_back(tube_segs);

	beg = end;

	++tube_id;
}

void read_occlu(string filepath, std::unordered_map<int, std::vector<cv::Point3i>>& orig_colli, std::unordered_map<int, std::set<int>>& occlu_graph)
{
	string filename = filepath + "/occlu.txt";
	std::ifstream infile(filename);
	std::string line;

	std::getline(infile, line);
	int number2;
	std::istringstream(line) >> number2;							// the number of all collisions


	for (int i = 0; i < number2; i++)
	{
		//std::cout << line << std::endl;
		std::getline(infile, line);
		int a, b, c, d, e, f;
		std::istringstream(line) >> a >> b >> c >> d >> e >> f;	// a: tube1, b: tube2, (c,d): the collision segment (beg->end)

		orig_colli[a].push_back(cv::Point3i(b, c, d));

		occlu_graph[a].insert(b);
	}
}