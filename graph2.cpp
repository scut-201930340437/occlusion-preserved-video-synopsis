#include "graph.h"
#include <algorithm>

void find_all_collide_segs(Segment* seg, std::vector<Segment*>& collide_segs)
{
	seg->is_traversed = true;
	collide_segs.push_back(seg);

	for (size_t i = 0; i < seg->colli_segs_.size(); i++)
	{
		if (!seg->colli_segs_[i]->is_traversed)
		{
			find_all_collide_segs(seg->colli_segs_[i], collide_segs);
		}
	}
}

/* Find the segment containing pose, and then split it and its neighbors*/
void find_and_split_segment_using_a(std::vector<Segment*>& heads, Segment* head, int pos)
{
	Segment* seg = head->next_;

	while (seg != NULL)
	{
		int a = seg->a_;
		int b = seg->b_;

		if (pos == a)
		{
			return;
		}

		if (pos > a && pos <= b)
		{
			std::vector<Segment*> collide_segs;
			reset_traverse(heads);
			find_all_collide_segs(seg, collide_segs);
			reset_traverse(heads);

			
			for (size_t i = 0; i < collide_segs.size(); i++)
			{
				Segment* seg1, * seg2;
				split_segment_using_a(collide_segs[i], pos, seg1, seg2);
				
				collide_segs[i]->divide_seg1_ = seg1;
				collide_segs[i]->divide_seg2_ = seg2;
				

				seg1->next_ = seg2;
				seg2->last_ = seg1;

				seg1->last_ = collide_segs[i]->last_;
				seg2->next_ = collide_segs[i]->next_;

				collide_segs[i]->last_->next_ = seg1;
				if (collide_segs[i]->next_ != NULL)
				{
					collide_segs[i]->next_->last_ = seg2;
				}
			}

			for (size_t i = 0; i < collide_segs.size(); i++)
			{
				Segment* tmp = collide_segs[i];
				for (size_t j = 0; j < tmp->colli_segs_.size(); j++)
				{
					Segment* tmp2 = tmp->colli_segs_[j];
					tmp->divide_seg1_->colli_segs_.push_back(tmp2->divide_seg1_);
					tmp->divide_seg2_->colli_segs_.push_back(tmp2->divide_seg2_);
				}
			}

			for (size_t i = 0; i < collide_segs.size(); i++)
			{
				delete collide_segs[i];
			}

			collide_segs.clear();

			break;
		}

		seg = seg->next_;
	}


	// The following is not right
	/*Segment* seg = head->next_;

	while (seg != NULL)
	{
		int a = seg->a_;
		int b = seg->b_;

		if (pos == a)
		{
			return;
		}

		if (pos > a && pos <= b)
		{
			Segment* seg1, * seg2;
			split_segment_using_a(seg, pos, seg1, seg2);

			seg1->next_ = seg2;
			seg2->last_ = seg1;

			seg1->last_ = seg->last_;
			seg2->next_ = seg->next_;

			seg->last_->next_ = seg1;
			if (seg->next_ != NULL)
			{
				seg->next_->last_ = seg2;
			}

			reset_traverse(heads);

			seg->is_traversed = true;
			for (size_t i = 0; i < seg->colli_segs_.size(); i++)
			{
				if (!seg->colli_segs_[i]->is_traversed)
				{
					seg->colli_segs_[i]->split_using_a(pos, seg1, seg2);
					delete seg->colli_segs_[i];
				}
			}

			delete seg;

			return;
		}

		seg = seg->next_;
	}*/
}

void find_and_split_segment_using_b(std::vector<Segment*>& heads, Segment* head, int pos)
{
	Segment* seg = head->next_;

	while (seg != NULL)
	{
		int a = seg->a_;
		int b = seg->b_;

		if (pos == b)
		{
			return;
		}

		if (pos >= a && pos < b)
		{
			std::vector<Segment*> collide_segs;
			reset_traverse(heads);
			find_all_collide_segs(seg, collide_segs);
			reset_traverse(heads);

			for (size_t i = 0; i < collide_segs.size(); i++)
			{
				Segment* seg1, * seg2;
				split_segment_using_b(collide_segs[i], pos, seg1, seg2);
				collide_segs[i]->divide_seg1_ = seg1;
				collide_segs[i]->divide_seg2_ = seg2;

				seg1->next_ = seg2;
				seg2->last_ = seg1;

				seg1->last_ = collide_segs[i]->last_;
				seg2->next_ = collide_segs[i]->next_;

				collide_segs[i]->last_->next_ = seg1;
				if (collide_segs[i]->next_ != NULL)
				{
					collide_segs[i]->next_->last_ = seg2;
				}
			}

			for (size_t i = 0; i < collide_segs.size(); i++)
			{
				Segment* tmp = collide_segs[i];
				for (size_t j = 0; j < tmp->colli_segs_.size(); j++)
				{
					Segment* tmp2 = tmp->colli_segs_[j];
					tmp->divide_seg1_->colli_segs_.push_back(tmp2->divide_seg1_);
					tmp->divide_seg2_->colli_segs_.push_back(tmp2->divide_seg2_);
				}
			}

			for (size_t i = 0; i < collide_segs.size(); i++)
			{
				delete collide_segs[i];
			}

			collide_segs.clear();
			break;
		}

		seg = seg->next_;
	}


	/*Segment* seg = head->next_;

	while (seg != NULL)
	{
		int a = seg->a_;
		int b = seg->b_;

		if (pos == b)
		{
			return;
		}

		if (pos >= a && pos < b)
		{
			Segment* seg1, * seg2;
			split_segment_using_b(seg, pos, seg1, seg2);

			seg1->next_ = seg2;
			seg2->last_ = seg1;

			seg1->last_ = seg->last_;
			seg2->next_ = seg->next_;

			seg->last_->next_ = seg1;
			if (seg->next_ != NULL)
			{
				seg->next_->last_ = seg2;
			}

			reset_traverse(heads);

			seg->is_traversed = true;
			for (size_t i = 0; i < seg->colli_segs_.size(); i++)
			{
				if (!seg->colli_segs_[i]->is_traversed)
				{
					seg->colli_segs_[i]->split_using_b(pos, seg1, seg2);
					delete seg->colli_segs_[i];
				}
			}

			delete seg;

			return;
		}

		seg = seg->next_;
	}*/
}

void ShowTubeSegmentsAndCollisions(Segment* head)
{
	Segment* next = head->next_;
	
	while (next != NULL)
	{
		int a = next->a_;
		int b = next->b_;

		for (int i = a; i <= b; i++)
		{

		}
	}
}

void find_circle(std::vector<Segment*> & heads, 
	Segment* begin, 
	std::vector<std::vector<Segment*>>& paths)
{
	std::stack<Segment*> Stack;
	
	Stack.push(begin);

	if (begin->next_ != NULL)
	{
		reset_traverse(heads);
		find_circle(begin, begin->next_, Stack, paths);
	}
	
	while (!Stack.empty())
	{
		Stack.pop();
	}
	Stack.push(begin);

	if (begin->last_ != NULL)
	{
		reset_traverse(heads);
		find_circle(begin, begin->last_, Stack, paths);
	}

	for (size_t i = 0; i < begin->colli_segs_.size(); i++)
	{
		while (!Stack.empty())
		{
			Stack.pop();
		}
		Stack.push(begin);

		reset_traverse(heads);
		find_circle(begin, begin->colli_segs_[i], Stack, paths);
	}
}

void find_circle(Segment* begin, Segment* current, 
	std::stack<Segment*>& stack,
	std::vector<std::vector<Segment*>> & paths)
{
	current->is_traversed = true;
	stack.push(current);

	//std::cout << "("<<current->tube_id_ << "," << current->a_ << "," << current->b_ << "),";
	
	if (current == begin)
	{
		std::vector<Segment*> path;
		while (!stack.empty())
		{
			path.push_back(stack.top());
			stack.pop();
		}

		//std::cout << "one path found" << std::endl;

		if(path.size() > 3)
			paths.push_back(path);

		for (int i = int(path.size()) - 1; i >=0; i--)
		{
			stack.push(path[i]);
		}
	}
	else
	{
		if (current->next_ != NULL && !current->next_->is_traversed && !current->is_next_break_)
		{
			find_circle(begin, current->next_, stack, paths);
		}

		if (current->last_ != NULL && !current->last_->is_traversed && !current->is_last_break_)
		{
			find_circle(begin, current->last_, stack, paths);
		}

		for (size_t i = 0; i < current->colli_segs_.size(); i++)
		{
			if (!current->colli_segs_[i]->is_traversed)
			{
				find_circle(begin, current->colli_segs_[i], stack, paths);
			}
		}
	}

	stack.pop();
	current->is_traversed = false;
}

void find_path(Segment* end, Segment* current,
	std::stack<Segment*>& stack,
	std::vector<std::vector<Segment*>>& paths)
{
	current->is_traversed = true;
	stack.push(current);

	if (current == end)
	{
		std::vector<Segment*> path;
		while (!stack.empty())
		{
			path.push_back(stack.top());
			stack.pop();
		}

		paths.push_back(path);

		for (int i = int(path.size()) - 1; i >= 0; i--)
		{
			stack.push(path[i]);
		}
	}
	else
	{
		if (current->next_ != NULL && !current->next_->is_traversed && !current->is_next_break_)
		{
			find_path(end, current->next_, stack, paths);
		}

		if (current->last_ != NULL && !current->last_->is_traversed && !current->is_last_break_)
		{
			find_path(end, current->last_, stack, paths);
		}

		for (size_t i = 0; i < current->colli_segs_.size(); i++)
		{
			if (!current->colli_segs_[i]->is_traversed)
			{
				find_path(end, current->colli_segs_[i], stack, paths);
			}
		}
	}

	stack.pop();
	current->is_traversed = false;
}

void dfs_all_connected_segments(Segment* node, std::vector<Segment*>& nodes)
{
	node->is_traversed = true;
	nodes.push_back(node);

	if (node->next_ != NULL && node->next_->is_traversed == false)
	{
		dfs_all_connected_segments(node->next_, nodes);
	}

	if (node->last_ != NULL && node->last_->is_traversed == false)
	{
		dfs_all_connected_segments(node->last_, nodes);
	}

	for (size_t i = 0; i < node->colli_segs_.size(); i++)
	{
		if (node->colli_segs_[i]->is_traversed == false)
		{
			dfs_all_connected_segments(node->colli_segs_[i], nodes);
		}
	}
}

void find_path_without_backtrack(Segment* end, Segment* current,
	std::stack<Segment*>& stack,
	std::vector<Segment*>& path)
{
	current->is_traversed = true;
	stack.push(current);

	if (current == end)
	{
		while (!stack.empty())
		{
			path.push_back(stack.top());
			stack.pop();
		}
	}
	else
	{
		if (path.empty() && current->next_ != NULL && !current->next_->is_traversed && !current->is_next_break_)
		{
			find_path_without_backtrack(end, current->next_, stack, path);
		}

		if (path.empty() && current->last_ != NULL && !current->last_->is_traversed && !current->is_last_break_)
		{
			find_path_without_backtrack(end, current->last_, stack, path);
		}

		for (size_t i = 0; i < current->colli_segs_.size(); i++)
		{
			if (path.empty() && !current->colli_segs_[i]->is_traversed && !current->colli_break_[i])
			{
				find_path_without_backtrack(end, current->colli_segs_[i], stack, path);
			}
		}
	}

	if (path.empty())
		stack.pop();
}

void find_path_without_backtrack_clockwise(Segment* end, Segment* current,
	std::stack<Segment*>& stack,
	std::vector<Segment*>& path)
{
	current->is_traversed = true;
	stack.push(current);

	if (current == end)
	{
		while (!stack.empty())
		{
			path.push_back(stack.top());
			stack.pop();
		}
	}
	else
	{
		for (size_t i = 0; i < current->colli_segs_.size(); i++)
		{
			if (path.empty() && !current->colli_segs_[i]->is_traversed && !current->colli_break_[i] 
				&& current->colli_segs_[i]->tube_id_ < current->tube_id_)
			{
				find_path_without_backtrack_clockwise(end, current->colli_segs_[i], stack, path);
			}
		}

		if (path.empty() && current->next_ != NULL && !current->next_->is_traversed && !current->is_next_break_)
		{
			find_path_without_backtrack_clockwise(end, current->next_, stack, path);
		}

		for (size_t i = 0; i < current->colli_segs_.size(); i++)
		{
			if (path.empty() && !current->colli_segs_[i]->is_traversed && !current->colli_break_[i]
				&& current->colli_segs_[i]->tube_id_ > current->tube_id_)
			{
				find_path_without_backtrack_clockwise(end, current->colli_segs_[i], stack, path);
			}
		}

		if (path.empty() && current->last_ != NULL && !current->last_->is_traversed && !current->is_last_break_)
		{
			find_path_without_backtrack_clockwise(end, current->last_, stack, path);
		}
	}

	if(path.empty())
		stack.pop();
}



void find_path_without_backtrack_anti_clockwise(Segment* end, Segment* current,
	std::stack<Segment*>& stack,
	std::vector<Segment*>& path)
{
	current->is_traversed = true;
	stack.push(current);

	if (current == end)
	{
		while (!stack.empty())
		{
			path.push_back(stack.top());
			stack.pop();
		}
	}
	else
	{
		if (path.empty() && current->last_ != NULL && !current->last_->is_traversed && !current->is_last_break_)
		{
			find_path_without_backtrack_anti_clockwise(end, current->last_, stack, path);
		}

		for (int i = int(current->colli_segs_.size()) - 1; i >= 0; i--)
		{
			if (path.empty() && !current->colli_segs_[i]->is_traversed && !current->colli_break_[i]
				&& current->colli_segs_[i]->tube_id_ > current->tube_id_)
			{
				find_path_without_backtrack_anti_clockwise(end, current->colli_segs_[i], stack, path);
			}
		}

		if (path.empty() && current->next_ != NULL && !current->next_->is_traversed && !current->is_next_break_)
		{
			find_path_without_backtrack_anti_clockwise(end, current->next_, stack, path);
		}

		for (int i = int(current->colli_segs_.size()) - 1; i >= 0; i--)
		{
			if (path.empty() && !current->colli_segs_[i]->is_traversed && !current->colli_break_[i]
				&& current->colli_segs_[i]->tube_id_ < current->tube_id_)
			{
				find_path_without_backtrack_anti_clockwise(end, current->colli_segs_[i], stack, path);
			}
		}
	}

	if (path.empty())
		stack.pop();
}



void break_links_in_path(std::vector<Segment*> path) {
	for (size_t i = 0; i < path.size() - 1; i++) {
		if (path[i]->tube_id_ == path[i + 1]->tube_id_) {
			if (path[i]->last_ == path[i + 1]) {
				path[i]->is_last_break_ = true;
				path[i + 1]->is_next_break_ = true;
			}
			else {
				path[i]->is_next_break_ = true;
				path[i + 1]->is_last_break_ = true;
			}
		}
		else {
			for (size_t j = 0; j < path[i]->colli_segs_.size(); j++) {
				if (path[i]->colli_segs_[j] == path[i + 1]) {
					path[i]->colli_break_[j] = true;
					break;
				}
			}
			for (size_t j = 0; j < path[i + 1]->colli_segs_.size(); j++) {
				if (path[i + 1]->colli_segs_[j] == path[i]) {
					path[i + 1]->colli_break_[j] = true;
					break;
				}
			}
		}
	}
}

bool find_circle_by_two_path(std::vector<Segment*> & heads, Segment* A, Segment* B, std::vector<Segment*>& path1, std::vector<Segment*>& path2)
{
	reset_last_next_break(heads);
	reset_traverse(heads);
	reset_colli_break(heads);

	std::stack<Segment*> tmp1;
	find_path_without_backtrack_clockwise(B, A, tmp1, path1);
	if (path1.empty())
		return false;

	break_links_in_path(path1);
	reset_traverse(heads);
	
	std::stack<Segment*> tmp2;
	find_path_without_backtrack_anti_clockwise(B, A, tmp2, path2);
	if (path2.empty())
		return false;

	break_links_in_path(path2);

	return true;
}

void find_rest_fish(std::vector<Segment*>& heads, std::vector<Segment*> fish, std::set<Segment*>& asic,int circulation_num)
{
	std::vector<Segment*> rest_fish;
	for (size_t i = 0; i < fish.size(); i++) {
		std::vector<Segment*> path1, path2;
		bool founded;
		

		if (fish[i]->next_ != NULL) 
		{
			/*bool isSameCir = false;
			for (std::set<int>::iterator iter1=fish[i]->circulation_id_.begin();iter1!=fish[i]->circulation_id_.end();++iter1)
			{
				if (fish[i]->next_->circulation_id_.count(*iter1)>0)
				{
					isSameCir = true;
					break;
				}
			}*/

			if (fish[i]->next_->circulation_id_.count(circulation_num) == 0)
			{
				path1.clear();
				path2.clear();
				founded = find_circle_by_two_path(heads, fish[i], fish[i]->next_, path1, path2);

				if (founded)
				{
					for (size_t j = 0; j < path1.size(); j++)
					{
						path1[j]->circulation_id_.insert(circulation_num);
						asic.insert(path1[j]);
					}

					for (size_t j = 0; j < path2.size(); j++)
					{
						path2[j]->circulation_id_.insert(circulation_num);
						asic.insert(path2[j]);
					}
					rest_fish.push_back(fish[i]->next_);
				}
			}
		}

		
		if (fish[i]->last_ != NULL) 
		{
			/*bool isSameCir = false;
			for (std::set<int>::iterator iter1 = fish[i]->circulation_id_.begin(); iter1 != fish[i]->circulation_id_.end(); ++iter1)
			{
				if (fish[i]->last_->circulation_id_.count(*iter1)>0)
				{
					isSameCir = true;
					break;
				}
			}*/

			if (fish[i]->last_->circulation_id_.count(circulation_num) == 0)
			{ 
				path1.clear();
				path2.clear();
				founded = find_circle_by_two_path(heads, fish[i], fish[i]->last_, path1, path2);

				if (founded)
				{
					for (size_t j = 0; j < path1.size(); j++)
					{
						path1[j]->circulation_id_.insert(circulation_num);
						asic.insert(path1[j]);
					}

					for (size_t j = 0; j < path2.size(); j++)
					{
						path2[j]->circulation_id_.insert(circulation_num);
						asic.insert(path2[j]);
					}
					rest_fish.push_back(fish[i]->last_);
				}
			}
		}

		//
		
		if (fish[i]->colli_segs_.size()) 
		{
			for (size_t k = 0; k < fish[i]->colli_segs_.size(); k++) 
			{
				/*bool isSameCir = false;
				for (std::set<int>::iterator iter1 = fish[i]->circulation_id_.begin(); iter1 != fish[i]->circulation_id_.end(); ++iter1)
				{
					if (fish[i]->colli_segs_[k]->circulation_id_.count(*iter1)>0)
					{
						isSameCir = true;
						break;
					}
				}*/

				if (fish[i]->colli_segs_[k]->circulation_id_.count(circulation_num) == 0) 
				{
					path1.clear();
					path2.clear();
					founded = find_circle_by_two_path(heads, fish[i], fish[i]->colli_segs_[k], path1, path2);

					if (founded)
					{
						for (size_t j = 0; j < path1.size(); j++)
						{
							path1[j]->circulation_id_.insert(circulation_num);
							asic.insert(path1[j]);
						}

						for (size_t j = 0; j < path2.size(); j++)
						{
							path2[j]->circulation_id_.insert(circulation_num);
							asic.insert(path2[j]);
						}
						rest_fish.push_back(fish[i]->colli_segs_[k]);
					}
				}
			}
		}
	}

	if (rest_fish.size() == 0) {
		return;
	}

	find_rest_fish(heads, rest_fish, asic, circulation_num);
}


void gather_all_segments_in_circulation_by_new_method(std::vector<Segment*>& heads, Segment* A, std::set<Segment*>& asic,int circulation_num)
{
	std::vector<Segment*> all_connected_segs;
	reset_traverse(heads);
	dfs_all_connected_segments(A, all_connected_segs);

	//std::cout << "*****all connected seg*******" << std::endl;
	/*for (size_t i = 0; i < all_connected_segs.size(); i++)
	{
		std::cout << all_connected_segs[i]->tube_id_ << " " << all_connected_segs[i]->a_ << " " << all_connected_segs[i]->b_ << std::endl;
	}*/

	for (size_t i = 0; i < all_connected_segs.size(); i++)
	{
		Segment* B = all_connected_segs[i];

		if (B == A)
			continue;

		//if (B->circulation_traversed)
		//{
		//	//std::cout << "**** B ******  " << B->tube_id_ << " " << B->a_ << " " << B->b_ << std::endl;
		//	continue;
		//}

		bool isSameCir = false;
		for (std::set<int>::iterator iter=A->circulation_id_.begin();iter != A->circulation_id_.end();++iter)
		{
			if (B->circulation_id_.count(*iter) > 0)
			{
				isSameCir = true;
				break;
			}
			
		}
		if (isSameCir)
		{
			continue;
		}


		/*if (B->circulation_id_.count(circulation_num) > 0)
		{
			continue;
		}*/
		

		std::vector<Segment*> path1, path2;

		bool founded = find_circle_by_two_path(heads, A, B, path1, path2);

		if (founded)
		{
			for (size_t j = 0; j < path1.size(); j++)
			{
				path1[j]->circulation_id_.insert(circulation_num);
				asic.insert(path1[j]);
			}

			for (size_t j = 0; j < path2.size(); j++)
			{
				path2[j]->circulation_id_.insert(circulation_num);
				asic.insert(path2[j]);
			}
		}
	}

	std::set<Segment*>::iterator it;
	std::vector<Segment*> fish;
	for (it = asic.begin(); it != asic.end(); it++) 
	{
		fish.push_back(*it);
	}
	find_rest_fish(heads, fish, asic, circulation_num);
}

void gather_all_segments_in_circulation(std::vector<Segment*> & heads, Segment* seg, std::set<Segment*>& asic)
{
	std::vector<std::vector<Segment*>> Paths;

	find_circle(heads, seg, Paths);

	for (size_t i = 0; i < Paths.size(); i++)
	{
		for (size_t j = 0; j < Paths[i].size(); j++)
		{
			asic.insert(Paths[i][j]);
		}
	}
}

/* Here can be simplified  
 We do not have to find all_circulated_segs 
 Instead, we can find all candidates on the whole graph
 Procedures behand can help distinguish segments actually cannot be divided */
void find_segments_can_be_divided(std::vector<Segment*>& ans,
	std::vector<Segment*>& all_circulated_segs,
	Segment* slf,
	int pos)
{
	for (std::vector<Segment*>::iterator iter = all_circulated_segs.begin();
		iter != all_circulated_segs.end();
		iter++)
	{
		Segment* obj = *iter;
		if (obj == slf)
			continue;
		if (pos >= obj->a_ && pos <= obj->b_)
		{
			ans.push_back(obj);
		}
	}
}

void find_segments_can_be_divided_whole_graph(std::set<Segment*>& ans,
	std::vector<Segment*> & heads,
	Segment* slf,
	int pos)
{
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;

		while (next != NULL)
		{
			/*if (next == slf)
			{
				next = next->next_;
				continue;
			}*/

			bool isSameCir = false;
			for (std::set<int>::iterator iter = slf->circulation_id_.begin(); iter != slf->circulation_id_.end(); ++iter)
			{
				if (next->circulation_id_.count(*iter) > 0)
				{
					isSameCir = true;
					break;
				}

			}
			if (isSameCir && pos >= next->a_ && pos <= next->b_)
			{
				ans.insert(next);
			}
			next = next->next_;
		}
	}
}

struct StruDivide
{
	int tube_provide_pos;
	int tube_divided;
	int pos;
	int a_or_b; // 0: a, 1;b
};

bool check_if_having_link_already(Segment* seg1, Segment* seg2)
{
	bool f1 = false;

	for (size_t i = 0; i < seg1->colli_segs_.size(); i++)
	{
		if (seg1->colli_segs_[i] == seg2)
		{
			f1 = true;
			break;
		}
	}

	bool f2 = false;

	for (size_t i = 0; i < seg2->colli_segs_.size(); i++)
	{
		if (seg2->colli_segs_[i] == seg1)
		{
			f2 = true;
			break;
		}
	}

	assert(f1 == f2);
	
	//if (f1 && f2)
	if(f1)
		return true;
	else
		return false;
}

bool check_if_two_segments_can_be_linked(std::vector<Segment*>& heads, 
	Segment* seg1, 
	Segment* seg2)
{
	std::set<Segment*> cands; // 所有可以被划分的seg
	
	find_segments_can_be_divided_whole_graph(cands, heads, seg1, seg1->a_);
	find_segments_can_be_divided_whole_graph(cands, heads, seg1, seg1->b_);
	
	find_segments_can_be_divided_whole_graph(cands, heads, seg2, seg2->a_);
	find_segments_can_be_divided_whole_graph(cands, heads, seg2, seg2->b_);

	/* Break links */
	reset_last_next_break(heads);

	for (std::set<Segment*>::iterator iter = cands.begin();iter != cands.end();iter++)
	{
		Segment* seg = *iter;
		seg->is_next_break_ = true;
		seg->is_last_break_ = true;

		seg->last_->is_next_break_ = true;
		if (seg->next_ != NULL)
		{
			seg->next_->is_last_break_ = true;
		}
	}

	/* Check left and right connection */
	reset_traverse(heads);
	std::vector<Segment*> seg1_lefts;
	seg1->collect_left_neighbors(seg1_lefts);

	reset_traverse(heads);
	std::vector<Segment*> seg1_rights;
	seg1->collect_right_neighbors(seg1_rights);

	reset_traverse(heads);
	std::vector<Segment*> seg2_lefts;
	seg2->collect_left_neighbors(seg2_lefts);

	reset_traverse(heads);
	std::vector<Segment*> seg2_rights;
	seg2->collect_right_neighbors(seg2_rights);


	bool pre_flag = false;
	bool post_flag = false;

	for (size_t k = 0; k < seg1_lefts.size(); k++)
	{
		for (size_t kk = 0; kk < seg2_lefts.size(); kk++)
		{
			reset_traverse(heads);

			std::stack<Segment*> Stack;
			std::vector<Segment*> Path;
			find_path_without_backtrack(seg2_lefts[kk], seg1_lefts[k], Stack, Path);

			if (!Path.empty())
			{
				pre_flag = true;
				break;
			}
		}
		if (pre_flag)
		{
			break;
		}
	}

	for (size_t k = 0; k < seg1_rights.size(); k++)
	{
		for (size_t kk = 0; kk < seg2_rights.size(); kk++)
		{
			reset_traverse(heads);
			std::stack<Segment*> Stack;
			std::vector<Segment*> Path;
			find_path_without_backtrack(seg2_rights[kk], seg1_rights[k], Stack, Path);

			
			if (!Path.empty())
			{
				post_flag = true;
				break;
			}
		}
		if (post_flag)
		{
			break;
		}
	}


	/*for (std::set<Segment*>::iterator iter = cands.begin(); iter != cands.end(); iter++)
	{
		Segment* seg = *iter;
		seg->is_next_break_ = false;
		seg->is_last_break_ = false;

		seg->last_->is_next_break_ = false;
		if (seg->next_ != NULL)
		{
			seg->next_->is_last_break_ = false;
		}
	}*/



	reset_last_next_break(heads);

	if (pre_flag && post_flag)
		return true;
	else
		return false;
}

void add_link_between_segments(Segment* seg1, Segment* seg2)
{
	bool f = false;
	for (size_t j = 0; j < seg1->colli_segs_.size(); j++)
	{
		if (seg1->colli_segs_[j] == seg2)
		{
			f = true;
			break;
		}
	}

	if (!f)
	{
		seg1->colli_segs_.push_back(seg2);
	}

	f = false;
	for (size_t j = 0; j < seg2->colli_segs_.size(); j++)
	{
		if (seg2->colli_segs_[j] == seg1)
		{
			f = true;
			break;
		}
	}

	if (!f)
	{
		seg2->colli_segs_.push_back(seg1);
	}
}



void remove_link_between_segments(Segment* seg1, Segment* seg2)
{
	std::vector<Segment*>::iterator iter = seg1->colli_segs_.begin();
	while (iter != seg1->colli_segs_.end())
	{
		if (*iter == seg2)
		{
			iter = seg1->colli_segs_.erase(iter);
			break;
		}
		else
		{
			++iter;
		}
	}

	iter = seg2->colli_segs_.begin();
	while (iter != seg2->colli_segs_.end())
	{
		if (*iter == seg1)
		{
			iter = seg2->colli_segs_.erase(iter);
			break;
		}
		else
		{
			++iter;
		}
	}
}

void add_links_after_add_key_link(std::vector<std::pair<Segment*, Segment*>>& cands, 
	std::vector<Segment*>& heads, std::vector<size_t> & idx)
{
	for (size_t i = 0; i < cands.size(); i++)
	{
		Segment* seg1 = cands[i].first;
		Segment* seg2 = cands[i].second;

		//bool flg = ;

		if (!check_if_having_link_already(seg1, seg2))
		{
			if (check_if_two_segments_can_be_linked(heads, seg1, seg2))
			{
				idx.push_back(i);
				add_link_between_segments(seg1, seg2);
				reset_colli_break(heads);
			}
		}
	}
}

bool check_left_and_right(std::vector<Segment*> & heads, Segment * seg1, Segment * seg2)
{
	reset_traverse(heads);
	std::vector<Segment*> seg1_lefts;
	seg1->collect_left_neighbors(seg1_lefts);

	reset_traverse(heads);
	std::vector<Segment*> seg1_rights;
	seg1->collect_right_neighbors(seg1_rights);

	reset_traverse(heads);
	std::vector<Segment*> seg2_lefts;
	seg2->collect_left_neighbors(seg2_lefts);

	reset_traverse(heads);
	std::vector<Segment*> seg2_rights;
	seg2->collect_right_neighbors(seg2_rights);


	bool pre_flag = false;
	bool post_flag = false;

	for (size_t k = 0; k < seg1_lefts.size(); k++)
	{
		for (size_t kk = 0; kk < seg2_lefts.size(); kk++)
		{
			reset_traverse(heads);

			std::stack<Segment*> Stack;
			std::vector<Segment*> Path;
			find_path_without_backtrack(seg2_lefts[kk], seg1_lefts[k], Stack, Path);

			if (!Path.empty())
			{
				pre_flag = true;
				break;
			}
		}
		if (pre_flag)
		{
			break;
		}
	}

	for (size_t k = 0; k < seg1_rights.size(); k++)
	{
		for (size_t kk = 0; kk < seg2_rights.size(); kk++)
		{
			reset_traverse(heads);
			std::stack<Segment*> Stack;
			std::vector<Segment*> Path;
			find_path_without_backtrack(seg2_rights[kk], seg1_rights[k], Stack, Path);


			if (!Path.empty())
			{
				post_flag = true;
				break;
			}
		}
		if (post_flag)
		{
			break;
		}
	}

	if (pre_flag && post_flag)
		return true;
	else
		return false;
}


void classify_links_after_dividing(std::vector<Segment*>& heads, std::vector<std::pair<Segment*, Segment*>>& links)
{
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;

		while (next != NULL)
		{
			if (next->circulation_id_.empty())
			{
				next = next->next_;
				continue;
			}

			std::set<Segment*> cands;

			find_segments_can_be_divided_whole_graph(cands, heads, next, next->a_);
			find_segments_can_be_divided_whole_graph(cands, heads, next, next->b_);

			reset_last_next_break(heads);

			/* Break links */
			for (std::set<Segment*>::iterator iter = cands.begin(); iter != cands.end(); iter++)
			{
				Segment* seg = *iter;
				seg->is_next_break_ = true;
				seg->is_last_break_ = true;

				seg->last_->is_next_break_ = true;
				if (seg->next_ != NULL)
				{
					seg->next_->is_last_break_ = true;
				}
			}

			for (std::set<Segment*>::iterator iter = cands.begin(); iter != cands.end(); iter++)
			{
				Segment* current = *iter;

				if (current == next)
				{
					continue;
				}

				/*if (check_left_and_right(heads, next, current))
				{
					bool kkkflag = false;
					for (size_t kkk = 0; kkk < links.size(); kkk++)
					{
						if ((next == links[kkk].first && current == links[kkk].second) ||
							(next == links[kkk].second && current == links[kkk].first))
						{
							kkkflag = true;
							break;
						}
					}
					if(!kkkflag)
						links.push_back(std::pair<Segment*, Segment*>(next, current));
				}*/
				if (!check_if_having_link_already(next, current))
				{
					bool kkkflag = false;
					for (size_t kkk = 0; kkk < links.size(); kkk++)
					{
						if ((next == links[kkk].first && current == links[kkk].second) ||
							(next == links[kkk].second && current == links[kkk].first))
						{
							kkkflag = true;
							break;
						}
					}
					if (!kkkflag)
						links.push_back(std::pair<Segment*, Segment*>(next, current));
				}
			}

			/*for (std::set<Segment*>::iterator iter = cands.begin(); iter != cands.end(); iter++)
			{
				Segment* seg = *iter;
				seg->is_next_break_ = false;
				seg->is_last_break_ = false;

				seg->last_->is_next_break_ = false;
				if (seg->next_ != NULL)
				{
					seg->next_->is_last_break_ = false;
				}
			}*/

			reset_last_next_break(heads);

			next = next->next_;
		}
	}
}


void find_divide_points( 
	Segment* current,
	std::vector<Segment*>& all_circulated_segs,
	std::vector<StruDivide>& divide_points,
	int a_or_b)
{

	std::vector<Segment*> divided_tubes;
	divided_tubes.push_back(current);

	int pos;

	if (a_or_b == 0)
	{
		pos = current->a_;
	}
	else
	{
		pos = current->b_;
	}

	find_segments_can_be_divided(divided_tubes, all_circulated_segs, current, pos);

	/*for (size_t i = 1; i < divided_tubes.size(); i++)
	{
		StruDivide sd;
		sd.tube_provide_pos = current->tube_id_;
		sd.tube_divided = divided_tubes[i]->tube_id_;
		sd.pos = pos;
		sd.a_or_b = a_or_b;
		divide_points.push_back(sd);
	}*/


	/* Break links */
	for (size_t i = 0; i < divided_tubes.size(); i++)
	{
		divided_tubes[i]->is_next_break_ = true;
		divided_tubes[i]->is_last_break_ = true;

		divided_tubes[i]->last_->is_next_break_ = true;
		if (divided_tubes[i]->next_ != NULL)
		{
			divided_tubes[i]->next_->is_last_break_ = true;
		}
	}

	/*reset_traverse(heads);
	std::vector<Segment*> current_lefts;
	current->collect_left_neighbors(current_lefts);

	reset_traverse(heads);
	std::vector<Segment*> current_rights;
	current->collect_right_neighbors(current_rights);*/


	for (size_t i = 1; i < divided_tubes.size(); i++)
	{
		//reset_traverse(heads);
		//std::vector<Segment*> objs_left;
		//divided_tubes[i]->collect_left_neighbors(objs_left);

		//reset_traverse(heads);
		//std::vector<Segment*> objs_right;
		//divided_tubes[i]->collect_right_neighbors(objs_right);

		//bool pre_flag = false;
		//bool post_flag = false;

		//for (size_t k = 0; k < current_lefts.size(); k++)
		//{
		//	for (size_t kk = 0; kk < objs_left.size(); kk++)
		//	{
		//		reset_traverse(heads);
		//
		//		std::stack<Segment*> Stack;
		//		std::vector<Segment*> Path;
		//		find_path_without_backtrack(objs_left[kk], current_lefts[k], Stack, Path);

		//		/*std::vector<std::vector<Segment*>> Paths;
		//		find_path(current_lefts[k], objs_left[kk], Stack, Paths);*/
		//		if (!Path.empty())
		//		{
		//			pre_flag = true;
		//			break;
		//		}
		//	}
		//	if (pre_flag)
		//	{
		//		break;
		//	}
		//}

		//for (size_t k = 0; k < current_rights.size(); k++)
		//{
		//	for (size_t kk = 0; kk < objs_right.size(); kk++)
		//	{
		//		reset_traverse(heads);
		//		std::stack<Segment*> Stack;
		//		std::vector<Segment*> Path;
		//		find_path_without_backtrack(objs_right[kk], current_rights[k], Stack, Path);

		//		/*std::vector<std::vector<Segment*>> Paths;
		//		find_path(current_rights[k], objs_right[kk], Stack, Paths);*/
		//		if (!Path.empty())
		//		{
		//			post_flag = true;
		//			break;
		//		}
		//	}
		//	if (post_flag)
		//	{
		//		break;
		//	}
		//}

		StruDivide sd;
		sd.tube_provide_pos = current->tube_id_;
		sd.tube_divided = divided_tubes[i]->tube_id_;
		sd.pos = pos;
		sd.a_or_b = a_or_b;

		/*if (pre_flag && post_flag)
		{
			divide_points.push_back(sd);
		}
		else
		{
			cannot_divide_points.push_back(sd);
		}*/
		divide_points.push_back(sd);
	}


	/*for (size_t i = 0; i < divided_tubes.size(); i++)
	{
		divided_tubes[i]->is_next_break_ = false;
		divided_tubes[i]->is_last_break_ = false;

		divided_tubes[i]->last_->is_next_break_ = false;
		if (divided_tubes[i]->next_ != NULL)
		{
			divided_tubes[i]->next_->is_last_break_ = false;
		}
	}*/

}

Segment* find_segment_at_pos(Segment* head, int pos)
{
	Segment* next = head->next_;
	while (next != NULL)
	{
		if (pos >= next->a_ && pos <= next->b_)
		{
			return next;
		}

		next = next->next_;
	}

	return NULL;
}


bool my_compare2(Segment* a, Segment* b)
{
	if (a->tube_id_ < b->tube_id_)
	{
		return true;
	}
	else if (a->tube_id_ > b->tube_id_)
	{
		return false;
	}
	else
	{
		if (a->a_ < b->a_)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}


bool check_consistency_stitch(std::vector<Segment*>& heads) // ends stitched together 
{
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;

		while (next != NULL)
		{
			Segment* nn = next->next_;

			if (nn != NULL)
			{
				if (next->bb_ + 1 != nn->aa_)
				{
					/*std::cout << "consistentcy stitch*********** " << next->tube_id_ << " " << next->a_ << " " << next->b_ << std::endl;
					std::cout << "consistentcy stitch    *********** " << next->tube_id_ << " " << next->bb_ << " " << nn->aa_ << std::endl;*/
					return false;
					/*std::cout << "check consistency wrong! ................................." << std::endl;
					std::cout << "(" << head->tube_id_ << ": " << head->a_ << "->" << head->b_ << "," << head->aa_ << "->" << head->bb_ << ")" << std::endl;
					std::cout << "(" << head->next_->tube_id_ << ": " << head->next_->a_ << "->" << head->next_->b_ << "," << head->next_->aa_ << "->" << head->next_->bb_ << ")" << std::endl;*/
				}
			}
			next = next->next_;
		}
	}
	return true;
}

struct saabb
{
	int aa;
	int bb;
};

void save_aa_bb(std::vector<saabb>& saved, std::vector<Segment*>& heads)
{
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;

		while (next != NULL)
		{
			saabb t;
			t.aa = next->aa_;
			t.bb = next->bb_;
			saved.push_back(t);

			next = next->next_;
		}
		
	}
}

void restore_aa_bb(std::vector<saabb>& saved, std::vector<Segment*>& heads)
{
	int idx = 0;
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;

		while (next != NULL)
		{
			next->aa_ = saved[idx].aa;
			next->bb_ = saved[idx].bb;
			idx++;
			next = next->next_;
		}
		
	}
}

void print_aabb(std::vector<Segment*>& heads)
{
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;

		while (next != NULL)
		{
			std::cout << next->tube_id_ << ",(" << next->a_ << "->" << next->b_ << ") (colli:";
			for (size_t j = 0; j < next->colli_segs_.size(); j++)
			{
				std::cout << next->colli_segs_[j]->tube_id_ << ",";
			}
			std::cout << ")" << std::endl;

			next = next->next_;
		}
	}
}

void compute_number_of_inconsistent_segments(std::vector<Segment*>& heads, std::vector<Segment*>& inconsis_segs)
{
	std::vector<saabb> saved;
	save_aa_bb(saved, heads);

	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;
	
		while (next != NULL)
		{
			

			int oaa = next->aa_;
			int obb = next->bb_;
			
			next->bb_ += 1000;

			reset_traverse(heads);
			
			dfs_change_aabb(next, oaa, obb);
			
			//std::cout << "change dfs comple****" << std::endl;
			reset_traverse(heads);

			//print_aabb(heads);
			
			if (!check_consistency_stitch(heads))
			{
				//std::cout << "check consis******  " <<" "<<next->tube_id_<<" "<<next->aa_<<" "<< next->last_->bb_ << std::endl;
				inconsis_segs.push_back(next);
			}
		
			restore_aa_bb(saved, heads);

			next = next->next_;
		}
	}
}

bool is_two_segments_linked(std::vector<Segment*>& heads, Segment* s1, Segment* s2)
{	
	reset_traverse(heads);
	std::vector<Segment*> collis;
	find_all_collide_segs(s1, collis);
	for (size_t i = 0; i < collis.size(); i++)
	{
		if (collis[i] == s2)
		{
			return true;
		}
	}
	return false;;
}

void resolve_circulation_problem(std::vector<Segment*> & heads)
{
	
	/* Find all dividing points for the whole graph */
	std::cout << "Find all dividing points for the whole graph" << std::endl;

	//std::vector<StruDivide> cannot_divide_points;
	// 
	// 搜索循环，计算冲突并解决 -> 重新搜索循环，计算冲突并解决（解决冲突时会加边，可能产生新的循环和新的冲突）
	while (true)
	{
		std::vector<StruDivide> divide_points;
		int circulation_num = 0;
		/*搜索出所有循环*/
		
		for (size_t i = 0; i < heads.size(); i++)
		{
			//std::cout << "Tube:" << i << std::endl;
			Segment* next = heads[i]->next_;

			while (next != NULL)
			{
				//std::cout << next->tube_id_ << "," << next->a_ << "," << next->b_ << std::endl;

				//if (next->circulation_traversed)
				//{
				//	/*std::cout << "cir_tra:true " << next->tube_id_ << " " << next->a_ << " " << next->b_ << std::endl;*/
				//	next = next->next_;
				//	continue;
				//}
				//next->circulation_traversed = true;

				std::set<Segment*> all_circulated_segs_set;

				//std::cout << "start gathering" << std::endl;
				//gather_all_segments_in_circulation(heads, next, all_circulated_segs);

				gather_all_segments_in_circulation_by_new_method(heads, next, all_circulated_segs_set, circulation_num);

				std::vector<Segment*> all_circulated_segs_vec;

				/*for (std::set<Segment*>::iterator iter = all_circulated_segs_set.begin(); iter != all_circulated_segs_set.end();
					iter++)
				{
					Segment* current = *iter;
					all_circulated_segs_vec.push_back(current);
				}*/
				all_circulated_segs_vec.assign(all_circulated_segs_set.begin(), all_circulated_segs_set.end());

				std::sort(all_circulated_segs_vec.begin(), all_circulated_segs_vec.end(), my_compare2); // 根据 tube id 和起始帧排序，

				/*for (std::set<Segment*>::iterator iter = all_circulated_segs.begin(); iter != all_circulated_segs.end();
					iter++)*/


				for (size_t iii = 0; iii < all_circulated_segs_vec.size(); iii++)
				{
					Segment* current = all_circulated_segs_vec[iii];
					current->circulation_id_.insert(circulation_num);
					//std::cout << "(" << current->tube_id_ << "," << current->a_ << "," << current->b_ << ")" << std::endl;

					reset_last_next_break(heads);
					//find_divide_points(heads, current, all_circulated_segs_vec, divide_points, cannot_divide_points, 0);
					find_divide_points(current, all_circulated_segs_vec, divide_points, 0);
					reset_last_next_break(heads);
					//find_divide_points(heads, current, all_circulated_segs_vec, divide_points, cannot_divide_points, 1);
					find_divide_points(current, all_circulated_segs_vec, divide_points, 1);

				}

				/*if (all_circulated_segs_vec.size() > 0)
				{
					std::cout << "******cir " << circulation_num << " **************" << std::endl;
					for (size_t iii = 0; iii < all_circulated_segs_vec.size(); iii++)
					{
						Segment* current = all_circulated_segs_vec[iii];
						std::cout << "(" << current->tube_id_ << "," << current->a_ << "," << current->b_ << ")" << " " << current->circulation_traversed << std::endl;
					}
				}*/

				//std::cout << "gathering ending" << std::endl;

				// 
				if (!all_circulated_segs_vec.empty())
				{
					// 下一个循环
					circulation_num++;
				}
				next = next->next_;
			}
		}
		
		/*计算所有冲突并解决*/
		/* Dividing */
		std::cout << "Dividing" << std::endl;
		std::cout << divide_points.size() << "***************Number of Dividing Points**************" << std::endl;
		for (size_t i = 0; i < divide_points.size(); i++)
		{
			//std::cout << i << std::endl;
			StruDivide sd = divide_points[i];
			if (sd.a_or_b == 0)
			{
				find_and_split_segment_using_a(heads, heads[sd.tube_divided], sd.pos);
			}
			else
			{
				find_and_split_segment_using_b(heads, heads[sd.tube_divided], sd.pos);
			}
		}

		// 统计further dividing之后，segment的数量
		int seg_num_after_resolve_conf = 0;
		for (size_t i = 0; i < heads.size(); i++)
		{
			Segment* head = heads[i];
			Segment* next = head->next_;

			while (next != NULL)
			{
				//std::cout << "(" << next->tube_id_ << ":" << next->a_ << "-" << next->b_;

				for (size_t j = 0; j < next->colli_segs_.size(); j++)
				{
					//std::cout << "," << next->colli_segs_[j]->tube_id_;
				}
				//std::cout << ")" << std::endl;
				next = next->next_;

				seg_num_after_resolve_conf++;
			}
			//std::cout << std::endl;
		}

		std::cout << "seg_num_after_resolve_conf: " << seg_num_after_resolve_conf << std::endl;

		reset_colli_break(heads);

		//std::vector<std::pair<Segment*, Segment*>> links, possible_links;
		std::vector<std::pair<Segment*, Segment*>> links;
		//classify_links_after_dividing(heads, links, possible_links);
		classify_links_after_dividing(heads, links);


		std::vector<Segment*> inconsis_segs;
		// 计算冲突数
		compute_number_of_inconsistent_segments(heads, inconsis_segs);

		std::cout << "inconsistent segments number: " << inconsis_segs.size() << std::endl;
		if (inconsis_segs.empty()) // -------------------------------------------------------------------------------
		{
			//std::cout << "inconsis_segs empty" << std::endl;
			break;
		}

		//
		/*std::cout << "**********inconsistent segments:**********" << std::endl;
		for (size_t t = 0; t < inconsis_segs.size(); t++)
		{
			std::cout << "tube id:" << inconsis_segs[t]->tube_id_ << " " << "beg->end:" << inconsis_segs[t]->aa_ << " " << inconsis_segs[t]->bb_ << std::endl;
		}*/
		//

		//
		/*std::cout << "*********possible links:***********" << std::endl;
		for (size_t t = 0; t < links.size(); t++)
		{
			std::cout << "tube1 id:" << links[t].first->tube_id_ << " " << "tube2 id:" << links[t].second->tube_id_ <<
				" beg->end:" << links[t].first->aa_ << " " << links[t].first->bb_ << " " << links[t].second->aa_ << " " << links[t].second->bb_ << std::endl;
		}*/
		//

		int min_inconsis = inconsis_segs.size();
		//std::vector<size_t> num_inconsis;
		int min_idx = -1;


		// 遍历，找到能使冲突下降最多的边 ----------------------------
		for (size_t i = 0; i < links.size(); i++)
		{
			std::cout << i << "/" << links.size() << std::endl;

			std::vector<size_t> idx;
			add_link_between_segments(links[i].first, links[i].second);
			reset_colli_break(heads);

			add_links_after_add_key_link(links, heads, idx);
			reset_colli_break(heads);

			inconsis_segs.clear();
			compute_number_of_inconsistent_segments(heads, inconsis_segs);

			// num_inconsis.push_back(inconsis_segs.size());
			if (inconsis_segs.size() < min_inconsis)
			{
				std::cout << "inconsis_segs decrease......" << std::endl;
				min_inconsis = inconsis_segs.size();
				std::cout << "inconsistent segs number:" << inconsis_segs.size() << std::endl;
				min_idx = int(i);
			}

			remove_link_between_segments(links[i].first, links[i].second);

			for (size_t j = 0; j < idx.size(); j++)
			{
				remove_link_between_segments(links[idx[j]].first, links[idx[j]].second);
			}
			reset_colli_break(heads);

			/*if (min_inconsis == 0 || inconsis_segs.size()<min_inconsis)
			{
				min_inconsis = inconsis_segs.size();
				break;
			}*/
			if (min_inconsis == 0)
			{
				break;
			}
		}

		// 将附加边实际添加进FOG -------------------
		std::vector<size_t> idx2;
		add_link_between_segments(links[min_idx].first, links[min_idx].second);
		reset_colli_break(heads);
		add_links_after_add_key_link(links, heads, idx2);
		reset_colli_break(heads);
		idx2.push_back(min_idx);


		//// 加上所有边
		//std::cout << "add all possible links....." << std::endl;
		//for (size_t i = 0; i < links.size(); i++)
		//{
		//	std::cout << i << "/" << links.size() << std::endl;

		//	std::vector<size_t> idx;
		//	add_link_between_segments(links[i].first, links[i].second);
		//	reset_colli_break(heads);

		//	add_links_after_add_key_link(links, heads, idx);
		//	reset_colli_break(heads);
		//}

		//// 计算冲突数
		//inconsis_segs.clear();
		//compute_number_of_inconsistent_segments(heads, inconsis_segs);
		//min_inconsis = int(inconsis_segs.size());
		
		if (min_inconsis == 0)
		{
			break;
		}

		// 清空circulation id，准备搜索新的循环
		for (size_t i = 0; i < heads.size(); ++i)
		{
			heads[i]->circulation_id_.clear();
			Segment* next = heads[i]->next_;
			while (next != NULL)
			{
				next->circulation_id_.clear();
				next = next->next_;
			}
		}
		

		// 把加入的边从links中删除
		/*int jj = 0;
		std::vector<std::pair<Segment*, Segment*>>::iterator iter = links.begin();
		while (iter != links.end())
		{
			bool fff = false;
			for (size_t uuu = 0; uuu < idx2.size(); uuu++)
			{
				if (idx2[uuu] == jj)
				{
					iter = links.erase(iter);
					fff = true;
					break;
				}
			}
			if (!fff)
			{
				++iter;
			}
			jj++;
		}*/
		


		///* Check how many possible links are not linked finally */
		//int stats_Num = 0;
		//for (size_t i = 0; i < links.size(); i++)
		//{
		//	Segment* s1 = links[i].first;
		//	Segment* s2 = links[i].second;

		//	if (!is_two_segments_linked(heads, s1, s2))
		//	{
		//		++stats_Num;
		//	}
		//}
		//std::cout << "====================================> Number of not linked:" << stats_Num << "/" << links.size() << std::endl;
		
	}
	
	//	if (sd.a_or_b == 0)
	//	{
	//		find_and_split_segment_using_a(heads, heads[sd.tube_divided], sd.pos);
	//	}
	//	else
	//	{
	//		find_and_split_segment_using_b(heads, heads[sd.tube_divided], sd.pos);
	//	}
	//}

	//std::cout << "**************links:******************" << std::endl;
	//for (size_t i = 0; i < links.size(); i++)
	//{
	//	// 将由于分割而消除的原有的links补回来
	//	add_link_between_segments(links[i].first, links[i].second);
	//	std::cout << links[i].first->tube_id_ << "," << links[i].second->tube_id_ << ",(" << links[i].first->a_ << "," << links[i].first->b_ << ")" << std::endl;
	//}
}



void resolve_circulation_problem_simply_adding_all_collisions(std::vector<Segment*>& heads)
{
	/* Find all dividing points for the whole graph */
	std::cout << "Find all dividing points for the whole graph" << std::endl;
	std::vector<StruDivide> divide_points;
	std::vector<StruDivide> cannot_divide_points;

	int circulation_num = 0;

	for (size_t i = 0; i < heads.size(); i++)
	{
		std::cout << "Tube:" << i << std::endl;
		Segment* next = heads[i]->next_;

		while (next != NULL)
		{
			//std::cout << next->tube_id_ << "," << next->a_ << "," << next->b_ << std::endl;
			if (next->circulation_traversed)
			{
				next = next->next_;
				continue;
			}
			next->circulation_traversed = true;

			std::set<Segment*> all_circulated_segs;

			//std::cout << "start gathering" << std::endl;
			//gather_all_segments_in_circulation(heads, next, all_circulated_segs);

			gather_all_segments_in_circulation_by_new_method(heads, next, all_circulated_segs, circulation_num);

			//std::cout << all_circulated_segs.size() << std::endl;

			/*for (std::set<Segment*>::iterator iter = all_circulated_segs.begin(); iter != all_circulated_segs.end();
				iter++)
			{
				Segment* current = *iter;
				std::cout << "(" << current->tube_id_ << "," << current->a_ << "," << current->b_ << ")" << std::endl;
			}*/



			std::vector<Segment*> all_circulated_segs_vec;

			for (std::set<Segment*>::iterator iter = all_circulated_segs.begin(); iter != all_circulated_segs.end();
				iter++)
			{
				Segment* current = *iter;
				all_circulated_segs_vec.push_back(current);
			}

			std::sort(all_circulated_segs_vec.begin(), all_circulated_segs_vec.end(), my_compare2);

			/*for (std::set<Segment*>::iterator iter = all_circulated_segs.begin(); iter != all_circulated_segs.end();
				iter++)*/

			for (size_t iii = 0; iii < all_circulated_segs_vec.size(); iii++)
			{
				Segment* current = all_circulated_segs_vec[iii];
				current->circulation_id_.insert(circulation_num);



				//std::cout << "(" << current->tube_id_ << "," << current->a_ << "," << current->b_ << ")" << std::endl;



				reset_last_next_break(heads);
				find_divide_points( current, all_circulated_segs_vec, divide_points, 0);
				reset_last_next_break(heads);
				find_divide_points( current, all_circulated_segs_vec, divide_points, 1);

				//current->circulation_traversed = true;
			}

			//std::cout << "gathering ending" << std::endl;

			if (!all_circulated_segs_vec.empty())
			{
				circulation_num++;
			}

			next = next->next_;
		}
	}

	/* Dividing */
	std::cout << "Dividing" << std::endl;

	/*std::vector<StruDivide> tttt;
	for (size_t i = 0; i < divide_points.size(); i++)
	{
		StruDivide sd = divide_points[i];
		if (sd.tube_provide_pos == 32)
		{
			tttt.push_back(sd);
		}
	}*/

	std::cout << divide_points.size() << "************* Number of Dividing Points ****************" << std::endl;
	for (size_t i = 0; i < divide_points.size(); i++)
	{
		std::cout << i << std::endl;

		StruDivide sd = divide_points[i];

		if (sd.a_or_b == 0)
		{
			find_and_split_segment_using_a(heads, heads[sd.tube_divided], sd.pos);
		}
		else
		{
			find_and_split_segment_using_b(heads, heads[sd.tube_divided], sd.pos);
		}
	}

	for (size_t i = 0; i < cannot_divide_points.size(); i++)
	{
		std::cout << i << std::endl;

		StruDivide sd = cannot_divide_points[i];



		if (sd.a_or_b == 0)
		{
			find_and_split_segment_using_a(heads, heads[sd.tube_divided], sd.pos);
		}
		else
		{
			find_and_split_segment_using_b(heads, heads[sd.tube_divided], sd.pos);
		}
	}

	reset_colli_break(heads);

	std::vector<std::pair<Segment*, Segment*>> links, possible_links;
	classify_links_after_dividing(heads, links);

	for (size_t i = 0; i < links.size(); i++)
	{
		add_link_between_segments(links[i].first, links[i].second);
		std::cout << links[i].first->tube_id_ << "," << links[i].second->tube_id_ << ",(" << links[i].first->a_ << "," << links[i].first->b_ << ")" << std::endl;

	}

	for (size_t i = 0; i < possible_links.size(); i++)
	{
		add_link_between_segments(possible_links[i].first, possible_links[i].second);
		//std::cout << possible_links[i].first->tube_id_ << "," << possible_links[i].second->tube_id_ << ",(" << possible_links[i].first->a_ << "," << possible_links[i].first->b_ << ")" << std::endl;
	}
}

bool my_compare(Segment* a, Segment* b)
{
	return (a->tube_id_ < b->tube_id_);
}

void divide_non_collide_segment_further(std::vector<Segment*>& heads, int min_div_len)
{
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;
		while (next != NULL)
		{
			if (next->colli_segs_.empty()) // this is a non-collision segment
			{
				int s = next->a_;
				int e = next->b_;
				int len = e - s + 1;
				int seg_num = len / min_div_len + 1; 
				int real_div_len = len / seg_num;
				
				for (int j = 0; j < seg_num - 1; j++)
				{
					Segment* seg1 = NULL, * seg2 = NULL;
					split_segment_using_b(next, next->a_ + real_div_len - 1, seg1, seg2);
					seg1->next_ = seg2;
					seg2->last_ = seg1;
					seg1->last_ = next->last_;
					seg2->next_ = next->next_;
					next->last_->next_ = seg1;
					if(next->next_!=NULL)
						next->next_->last_ = seg2;
					delete next;
					next = seg2;
				}
			}
			next = next->next_;
		}
	}
}

void build_full_collision_graph(std::vector<Segment*>& heads, int min_div_len)
{
	/* Find roots */
	std::vector<Segment*> roots;
	reset_traverse(heads);

	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* next = heads[i]->next_;
		while (next != NULL)
		{
			if (!next->is_traversed)
			{
				roots.push_back(next);
				dfs_traverse(next);
			}

			next = next->next_;
		}
	}

	/* Process each group */
	for (size_t i = 0; i < roots.size(); i++)
	{
		reset_traverse(heads);
		std::vector<Segment*> segs;
		dfs_collect_graph_segments(roots[i], segs);

		/* Get tubes */
		std::set<int> tube_ids;
		std::vector<int> tube_ids_vec;
		for (size_t j = 0; j < segs.size(); j++)
		{
			tube_ids.insert(segs[j]->tube_id_);
		}
		for (std::set<int>::iterator iter2 = tube_ids.begin(); iter2 != tube_ids.end(); iter2++)
		{
			tube_ids_vec.push_back(*iter2);
		}

		/* Get first and last frames, and all boundaries of segs */
		int max_frame = -1000000;
		int min_frame = 1000000;
		std::set<int> bdy;
		for (size_t j = 0; j < segs.size(); j++)
		{
			if (segs[j]->a_ < min_frame)
				min_frame = segs[j]->a_;
			if (segs[j]->b_ > max_frame)
				max_frame = segs[j]->b_;
			bdy.insert(segs[j]->b_);
		}
		
		/* Insert regular split points */
		int len = max_frame - min_frame + 1;
		int seg_num = len / min_div_len + 1;
		int real_div_len = len / seg_num;

		for (int j = 0; j < seg_num - 1; j++)
		{
			bdy.insert(min_frame + (j+1) * real_div_len - 1);
		}

		/* Divide */
		for (std::set<int>::iterator iter = bdy.begin(); iter != bdy.end(); iter++)
		{
			for (std::set<int>::iterator iter2 = tube_ids.begin(); iter2 != tube_ids.end(); iter2++)
			{
				find_and_split_segment_using_b(heads, heads[*iter2], *iter);
			}
		}
	
		/* Add collision relationships */
		for (size_t j = 0; j < tube_ids_vec.size(); j++)
		{
			Segment* tmp = heads[tube_ids_vec[j]]->next_;
			while (tmp != NULL)
			{
				for (size_t k = 0; k < tube_ids_vec.size(); k++)
				{
					if (k == j)
					{
						continue;
					}
					Segment* tmp2 = heads[tube_ids_vec[k]]->next_;
					while (tmp2 != NULL)
					{
						if (tmp->a_ == tmp2->a_ && tmp->b_ == tmp2->b_)
						{
							bool not_add = false;
							for (size_t l = 0; l < tmp->colli_segs_.size(); l++)
							{
								if (tmp->colli_segs_[l] == tmp2)
								{
									not_add = true;
									break;
								}
							}

							if (!not_add)
							{
								tmp->colli_segs_.push_back(tmp2);
							}
						}
						tmp2 = tmp2->next_;
					}

				}
				tmp = tmp->next_;
			}
		}
	}
}

cv::Mat create_segment_graph2(std::string filename,
	std::vector<Segment*>& heads,
	std::vector<std::vector<Segment*>>& all_segs,
	std::vector<std::pair<EdgeNode*, EdgeNode*>>& all_edges,
	int seg_len,
	int& srclen,
	int algo_type,
	bool non_colliding_dividing_component,
	std::unordered_map<int, std::vector<cv::Point3i>>&orig_colli)
{
	std::ifstream infile(filename);
	std::string line;

	std::getline(infile, line);
	int number1;			// number of tubes
	int max_frame;			// number of frames of src
	std::istringstream(line) >> number1;
	
	std::getline(infile, line);
	std::istringstream(line) >> max_frame;

	srclen = max_frame;

	std::vector<cv::Point2i> tube_ab;			// tubes' beg and end

	double average_tube_len = 0;

	for (int i = 0; i < number1; i++)
	{
		std::getline(infile, line);
		int no, a, b;
		std::istringstream(line) >> no >> a >> b;
		tube_ab.push_back(cv::Point2i(a, b));
		average_tube_len += (b - a + 1);
	}

	average_tube_len /= number1;
	std::cout << "Average tube length: " << average_tube_len << std::endl;

	std::vector<std::vector<cv::Point3i>> tube_colli_pre;		// collision positions between tubes
	std::vector<std::vector<cv::Point3i>> tube_colli;
	tube_colli_pre.resize(number1);								// store each tube's collision points
	std::getline(infile, line);
	int number2;
	std::istringstream(line) >> number2;							// the number of all collisions

	
	for (int i = 0; i < number2; i++)
	{
		//std::cout << line << std::endl;
		std::getline(infile, line);
		int a, b, c, d, e, f;
		std::istringstream(line) >> a >> b >> c >> d >> e >> f;	// a: tube1, b: tube2, (c,d): the collision segment (beg->end)
		tube_colli_pre[a - 1].push_back(cv::Point3i(b - 1, c, d));

		orig_colli[a].push_back(cv::Point3i(b, c, d));
		//ret.at<uchar>(a - 1, b - 1) = 1;
		//tube_colli[d - 1].push_back(cv::Point3i(a - 1, b, c));
	}

	tube_colli.resize(number1);
	for (size_t i = 0; i < tube_colli_pre.size(); i++)
	{
		if (tube_colli_pre[i].empty())
		{
			continue;
		}

		for (size_t j = 0; j < tube_colli_pre[i].size(); j++)
		{
			int id = tube_colli_pre[i][j].x;
			int aa = tube_colli_pre[i][j].y;
			int bb = tube_colli_pre[i][j].z;

			bool flag = false;

			for (size_t k = 0; k < tube_colli[i].size(); k++)
			{
				if (tube_colli[i][k].x == id)     // combine two segments of the same tube (把跟某个管（例如0号管）碰撞的所有段连起来)
				{
					if (aa < tube_colli[i][k].y)   // select smaller starting index
					{
						tube_colli[i][k].y = aa;
					}
					if (bb > tube_colli[i][k].z)   // select bigger ending index
					{
						tube_colli[i][k].z = bb;
					}
					flag = true;
					break;
				}
			}

			if (!flag)
			{
				tube_colli[i].push_back(cv::Point3i(id, aa, bb));
			}
		}
	}

	

	for (int i = 0; i < number1; i++)
	{
		Segment* head = new Segment(-1, -1, i, -1); // 伪head
		Segment* seg = new Segment(tube_ab[i].x, tube_ab[i].y, i, tube_ab[i].y - tube_ab[i].x + 1);	// 整个管
		head->next_ = seg;
		seg->last_ = head;
		heads.push_back(head);
	}
	
	/* Build collision graph */
	for (int i = 0; i < number1; i++)
	{
		std::cout << i << std::endl;
		for (size_t j = 0; j < tube_colli[i].size(); j++)
		{
			int id = tube_colli[i][j].x;
			int a = tube_colli[i][j].y;
			int b = tube_colli[i][j].z;

			// Find a segment in tube i that contains a
			// divide this segment into two segments, replacing it with the two segments
			find_and_split_segment_using_a(heads, heads[i], a);
			find_and_split_segment_using_b(heads, heads[i], b);

			// Find a segment in tube id that contains a
			// divide this segment into two segments, replacing it with the two segments
			find_and_split_segment_using_a(heads, heads[id], a);
			find_and_split_segment_using_b(heads, heads[id], b);

			Segment* head = heads[i];
			Segment* next = head->next_;

			// 多个碰撞重叠在一个物体上，在重叠边界处对遮挡段进行进一步划分
			while (next != NULL)
			{
				if (next->a_ >= a && next->b_ <= b)
				{
					find_and_split_segment_using_a(heads, heads[id], next->a_);
					find_and_split_segment_using_b(heads, heads[id], next->b_);
				}
				next = next->next_;
			}

			head = heads[id];
			next = head->next_;

			while (next != NULL)
			{
				if (next->a_ >= a && next->b_ <= b)
				{
					find_and_split_segment_using_a(heads, heads[i], next->a_);
					find_and_split_segment_using_b(heads, heads[i], next->b_);
				}
				next = next->next_;
			}
			

			std::vector<Segment*> s1, s2;
			head = heads[i];
			next = head->next_;

			while (next != NULL)
			{
				if (next->a_ >= a && next->b_ <= b)
				{
					s1.push_back(next);
				}
				next = next->next_;
			}

			head = heads[id];
			next = head->next_;

			while (next != NULL)
			{
				if (next->a_ >= a && next->b_ <= b)
				{
					s2.push_back(next);
				}
				next = next->next_;
			}

			assert(s1.size() == s2.size());

			for (size_t k = 0; k < s1.size(); k++)
			{
				bool flag = false;
				for (size_t kk = 0; kk < s1[k]->colli_segs_.size(); kk++)
				{
					if (s2[k] == s1[k]->colli_segs_[kk])
					{
						flag = true;
						break;
					}
				}
				if(!flag)
					s1[k]->colli_segs_.push_back(s2[k]);
				
				flag = false;
				for (size_t kk = 0; kk < s2[k]->colli_segs_.size(); kk++)
				{
					if (s1[k] == s2[k]->colli_segs_[kk])
					{
						flag = true;
						break;
					}
				}
				
				if(!flag)
					s2[k]->colli_segs_.push_back(s1[k]);
			}

			
		}
	}
	
	/* Sort collision segments 按照tube id 排序 */
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;

		while (next != NULL)
		{
			if(!next->colli_segs_.empty())
				std::sort(next->colli_segs_.begin(), next->colli_segs_.end(), my_compare);
			next = next->next_;
		}
	}

	int seg_num_before_resolve_conf = 0;

	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;

		while (next != NULL)
		{
			//std::cout << "(" << next->tube_id_ << ":" << next->a_ << "-" << next->b_;

			for (size_t j = 0; j < next->colli_segs_.size(); j++)
			{
				//std::cout << "," << next->colli_segs_[j]->tube_id_;
			}
			//std::cout << ")" << std::endl;
			next = next->next_;

			seg_num_before_resolve_conf++;
		}
		//std::cout << std::endl;
	}

	std::cout << "seg_num_before_resolve_conf: " << seg_num_before_resolve_conf << std::endl;

	/*std::vector<Segment*> path1, path2;
	bool r1 = find_circle_by_two_path(heads, heads[0]->next_, heads[0]->next_->next_, path1, path2);*/

	/*std::vector<Segment*> all_connected_segs;
	dfs_all_connected_segments(heads[0]->next_->next_, all_connected_segs);*/

	/*std::set<Segment*> asic;
	gather_all_segments_in_circulation_by_new_method(heads, heads[0]->next_->next_, asic);*/
		
	/*for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* node = heads[i];
		int cnt = 0;
		while (node != NULL)
		{
			cnt++;
			node = node->next_;
		}

		std::cout << cnt << ","<< tube_colli[i].size()<< std::endl;
	}*/


	if (algo_type == 0) // our method
	{
		/* Graph 1: consider circuit */
		/* Resolve circulation problem */
		
		std::cout << "Resolve circulation problem" << std::endl;
		resolve_circulation_problem(heads);

		/* Divide non-collide segments further */
		if(non_colliding_dividing_component)
			divide_non_collide_segment_further(heads, seg_len);
	}
	else if(algo_type == 1) 
	{
		/* Add full collisions in the collision region */
		resolve_circulation_problem_simply_adding_all_collisions(heads);
		
		if(non_colliding_dividing_component)
			divide_non_collide_segment_further(heads, seg_len);
	}
	else if (algo_type == 2)
	{
		/* Nie, Li */
		/* Or Graph 2: do not consider circuit */
		build_full_collision_graph(heads, seg_len);
	}

	/* all segs */
	all_segs.clear();
	int segi = 0;
	for (size_t i = 0; i < heads.size(); i++)
	{
		std::vector<Segment*> tmp;
		Segment* head = heads[i];
		head = head->next_;
		while (head != NULL)
		{
			head->seg_index_ = segi++; // segment index among all the segments
			tmp.push_back(head);
			head = head->next_;
		}

		all_segs.push_back(tmp);
	}

	std::cout << "The number of All segments:" << segi << std::endl;

	/* Used for computing collision */
	cv::Mat ret = cv::Mat(segi, segi, CV_8UC1, cv::Scalar(0)); // A matrix with element indicating whether two segments collide or not

	/* Fill the segment collision indication matrix */
	
	for (size_t i = 0; i < all_segs.size(); i++)
	{
		for (size_t j = 0; j < all_segs[i].size(); j++)
		{
			Segment* s = all_segs[i][j];
			int s_seg_i = s->seg_index_;
			for (size_t k = 0; k < s->colli_segs_.size(); k++)
			{
				int c_seg_i = s->colli_segs_[k]->seg_index_;
				ret.at<uchar>(s_seg_i, c_seg_i) = 1;
			}
		}
	}

	return ret;
}
