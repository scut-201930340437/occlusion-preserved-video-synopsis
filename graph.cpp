#include "graph.h"
#include <queue>


void get_front_segment(Segment* root, Segment*& front)
{
	root->is_traversed = true;
	if (front == NULL)
	{
		front = root;
	}
	else
	{
		if (root->a_ < front->a_)
		{
			front = root;
		}
	}

	EdgeNode* edge = root->firstedge;

	while (edge != NULL)
	{
		if (edge->is_retained && edge->seg->is_traversed == false && !edge->is_break)
			get_front_segment(edge->seg, front);

		edge = edge->next;
	}
}

void get_front_segment2(Segment* root, Segment*& front)
{
	root->is_traversed = true;

	if (root->last_ != NULL) // head 
	{
		if (front == NULL)
		{
			front = root;
		}
		else
		{
			if (root->a_ < front->a_)
			{
				front = root;
			}
		}
	}

	if (root->last_ != NULL && !root->last_->is_traversed)
	{
		get_front_segment2(root->last_, front);
	}

	if (root->next_ != NULL && !root->next_->is_traversed)
	{
		get_front_segment2(root->next_, front);
	}


	for (size_t i = 0; i < root->colli_segs_.size(); i++)
	{
		if (!root->colli_segs_[i]->is_traversed)
		{
			get_front_segment2(root->colli_segs_[i], front);
		}
	}
}

void get_back_segment(Segment* root, Segment*& back)
{
	root->is_traversed = true;
	if (back == NULL)
	{
		back = root;
	}
	else
	{
		if (root->b_ > back->b_)
		{
			back = root;
		}
	}

	EdgeNode* edge = root->firstedge;

	while (edge != NULL)
	{
		if (edge->is_retained && edge->seg->is_traversed == false && !edge->is_break)
			get_back_segment(edge->seg, back);

		edge = edge->next;
	}
}

void get_back_segment2(Segment* root, Segment*& back)
{
	root->is_traversed = true;

	
	if (back == NULL)
	{
		back = root;
	}
	else
	{
		if (root->b_ > back->b_)
		{
			back = root;
		}
	}


	if (root->last_ != NULL && !root->last_->is_traversed)
	{
		get_back_segment2(root->last_, back);
	}

	if (root->next_ != NULL && !root->next_->is_traversed)
	{
		get_back_segment2(root->next_, back);
	}


	for (size_t i = 0; i < root->colli_segs_.size(); i++)
	{
		if (!root->colli_segs_[i]->is_traversed)
		{
			get_back_segment2(root->colli_segs_[i], back);
		}
	}
}

void bfs_traverseMST_and_show(Context & context, Segment* seg, bool show, 
	std::vector<std::pair<Segment*,Segment*>> & influ_pairs)
{
	seg->is_traversed = true;
	if (show)
		std::cout << *seg;

	int tid = seg->tube_id_;
	int sid = 0;
	
	for (size_t i = 0; i < context.all_segs[tid].size(); i++)
	{
		if (context.all_segs[tid][i] == seg)
		{
			sid = int(i);
			break;
		}
	}

	for (int j = sid - 1; j >= 0; j--)
	{
		if (show)
			std::cout << *(context.all_segs[tid][j]);
		context.all_segs[tid][j]->is_traversed = true;

		influ_pairs.push_back(std::pair<Segment*, Segment*>(
			context.all_segs[tid][j+1], context.all_segs[tid][j]));
	}

	for (int j = sid + 1; j < int(context.all_segs[tid].size()); j++)
	{
		if (show)
			std::cout << *(context.all_segs[tid][j]);
		context.all_segs[tid][j]->is_traversed = true;

		influ_pairs.push_back(std::pair<Segment*, Segment*>(
			context.all_segs[tid][j - 1], context.all_segs[tid][j]));
	}

	for (size_t j = 0; j < context.all_segs[tid].size(); j++)
	{
		EdgeNode* edge = context.all_segs[tid][j]->firstedge;

		while (edge != NULL)
		{
			if (edge->is_retained && !edge->is_break 
				&& !edge->seg->is_traversed
				&& edge->seg->tube_id_ != context.all_segs[tid][j]->tube_id_)
			{
				influ_pairs.push_back(std::pair<Segment*, Segment*>(
					context.all_segs[tid][j], edge->seg));

				bfs_traverseMST_and_show(context, edge->seg, show, influ_pairs);
			}
			edge = edge->next;
		}
	}
}

void dfs_traverse(Segment* seg)
{
	seg->is_traversed = true;

	if (seg->next_ != NULL && !seg->next_->is_traversed)
	{
		dfs_traverse(seg->next_);
	}
	if (seg->last_ != NULL && !seg->last_->is_traversed)
	{
		dfs_traverse(seg->last_);
	}
	for (size_t i = 0; i < seg->colli_segs_.size(); i++)
	{
		if (!seg->colli_segs_[i]->is_traversed)
		{
			dfs_traverse(seg->colli_segs_[i]);
		}
	}
}


void dfs_collect_graph_segments(Segment* seg, std::vector<Segment*> & graph_segs)
{
	seg->is_traversed = true;
	graph_segs.push_back(seg);

	if (seg->next_ != NULL && !seg->next_->is_traversed)
	{
		dfs_collect_graph_segments(seg->next_,graph_segs);
	}
	if (seg->last_ != NULL && !seg->last_->is_traversed && seg->last_->last_ != NULL)
	{
		dfs_collect_graph_segments(seg->last_, graph_segs);
	}
	for (size_t i = 0; i < seg->colli_segs_.size(); i++)
	{
		if (!seg->colli_segs_[i]->is_traversed)
		{
			dfs_collect_graph_segments(seg->colli_segs_[i], graph_segs);
		}
	}
}


void extract_segment_from(int beg, int end, int tube_id, int seg_len, std::vector<Segment*>& segs)
{
	double num = double(end - beg + 1) / double(seg_len);
	if (num < 1)
	{
		num = 1;
	}

	int N = int(num + 0.5);

	int L = int(double(end - beg + 1) / double(N) + 0.5);

	std::vector<int> tmp;

	for (int i = 0; i < N; i++)
	{
		tmp.push_back(beg + i * L);
	}

	tmp.push_back(end+1);

	for (int i = 0; i < N; i++)
	{
		segs.push_back(new Segment(tmp[i], tmp[i + 1]-1, tube_id, tmp[i + 1] - tmp[i]));
	}
}
//
//void extract_segment_from(int beg, int end, int tube_id, int seg_len, std::vector<Segment*>& segs)
//{
//	int a = beg;
//	while (a + seg_len - 1 <= end)
//	{
//		int b = a + seg_len - 1;
//
//		Segment* seg = new Segment(a, b, tube_id, seg_len);
//
//		segs.push_back(seg);
//
//		a += seg_len;
//	}
//
//	if (a <= end)
//	{
//		Segment* seg = new Segment(a, end, tube_id, end - a + 1);
//		segs.push_back(seg);
//	}
//}

/* Get all segments of graph A */
void dfs_connected_graph_ofA(Segment* A, std::vector<Segment*>& all)
{
	A->is_traversed = true;
	all.push_back(A);

	EdgeNode* edge = A->firstedge;
	while (edge != NULL)
	{
		if (edge->is_retained && edge->seg->is_traversed == false && !edge->is_break)
		{
			dfs_connected_graph_ofA(edge->seg, all);
		}

		edge = edge->next;
	}
}

void dfs_jump_chain(Segment* A, std::vector<std::pair<int, int>>& chain)
{
	A->is_traversed = true;
	
	EdgeNode* edge = A->firstedge;

	while (edge != NULL)
	{
		if (edge->is_retained && edge->seg->is_traversed == false && !edge->is_break)
		{
			if (edge->seg->tube_id_ != A->tube_id_)
			{
				chain.push_back(std::pair<int, int>(A->tube_id_, edge->seg->tube_id_));
			}

			dfs_jump_chain(edge->seg, chain);
		}

		edge = edge->next;
	}
}

struct tmpNode
{
	int id;
	std::vector<tmpNode*> next;
	bool is_traverse;
};

void reset_tmpNode_traverse(std::vector<tmpNode*> & all_nodes)
{
	for (size_t i = 0; i < all_nodes.size(); i++)
	{
		all_nodes[i]->is_traverse = false;
	}
}

bool exist_path(tmpNode* node1, tmpNode* node2)
{
	node1->is_traverse = true;

	if (node1->id == node2->id)
	{
		return true;
	}
	else
	{
		for (size_t i = 0; i < node1->next.size(); i++)
		{
			if (node1->next[i]->is_traverse != true && exist_path(node1->next[i], node2))
			{
				return true;
			}
		}
	}
	return false;
}

bool check_circulation(std::vector<std::pair<int, int>>& chain)
{
	if (chain.empty())
	{
		return false;
	}

	
	std::vector<tmpNode*> all_nodes;

	for (size_t i = 0; i < chain.size(); i++)
	{
		bool exist = false;
		for (size_t j = 0; j < all_nodes.size(); j++)
		{
			if (all_nodes[j]->id == chain[i].first)
			{
				exist = true;
				break;
			}
		}
		if (!exist)
		{
			tmpNode* node = new tmpNode;
			node->id = chain[i].first;
			all_nodes.push_back(node);
		}

		exist = false;
		for (size_t j = 0; j < all_nodes.size(); j++)
		{
			if (all_nodes[j]->id == chain[i].second)
			{
				exist = true;
				break;
			}
		}
		if (!exist)
		{
			tmpNode* node = new tmpNode;
			node->id = chain[i].second;
			all_nodes.push_back(node);
		}
	}

	for (size_t i = 0; i < chain.size(); i++)
	{
		tmpNode* node1 = NULL, * node2 = NULL;
		for (size_t j = 0; j < all_nodes.size(); j++)
		{
			if (all_nodes[j]->id == chain[i].first)
			{
				node1 = all_nodes[j];
			}
			if (all_nodes[j]->id == chain[i].second)
			{
				node2 = all_nodes[j];
			}
		}

		node1->next.push_back(node2);
	}

	bool ret = false;
	for (size_t i = 0; i < all_nodes.size(); i++)
	{
		tmpNode* node = all_nodes[i];

		for (size_t j = 0; j < node->next.size(); j++)
		{
			reset_tmpNode_traverse(all_nodes);
			ret = exist_path(node->next[j], node);

			if (ret)
			{
				break;
			}
		}
		if (ret)
			break;
	}

	for (size_t i = 0; i < all_nodes.size(); i++)
	{
		delete all_nodes[i];
	}

	return ret;
}

bool check_connection(Segment* A, Segment* B)
{
	std::vector<Segment*> all;

	dfs_connected_graph_ofA(A, all);

	bool flag = false;

	for (size_t i = 0; i < all.size(); i++)
	{
		if (all[i] == B)
		{
			flag = true;
			break;
		}
	}

	return flag;
}

void reset_traverse(std::vector<std::vector<Segment*>>& all_segs)
{
	for (size_t i = 0; i < all_segs.size(); i++)
	{
		for (size_t j = 0; j < all_segs[i].size(); j++)
		{
			all_segs[i][j]->is_traversed = false;
		}
	}
}

void reset_traverse(std::vector<Segment*>& heads)
{
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		head->is_traversed = false;
		Segment* next = head->next_;
		while (next != NULL)
		{
			next->is_traversed = false;
			next = next->next_;
		}
	}
}

void reset_last_next_break(std::vector<Segment*>& heads)
{
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;
		while (next != NULL)
		{
			next->is_next_break_ = false;
			next->is_last_break_ = false;
			next = next->next_;
		}
	}
}

void reset_colli_break(std::vector<Segment*>& heads)
{
	for (size_t i = 0; i < heads.size(); i++)
	{
		Segment* head = heads[i];
		Segment* next = head->next_;
		while (next != NULL)
		{
			next->colli_break_.clear();
			next->colli_break_.resize(next->colli_segs_.size());
			for (size_t i = 0; i < next->colli_break_.size(); i++) {
				next->colli_break_[i] = false;
			}
			next = next->next_;
		}
	}
}


//void dfs_change_aabb(Segment* seg, int oaa, int obb)
//{
//	seg->is_traversed = true;
//
//	EdgeNode* edge = seg->firstedge;
//
//	while (edge != NULL)
//	{
//		if (edge->is_retained && !edge->seg->is_traversed && !edge->is_break)
//		{
//			Segment* adjseg = edge->seg;
//			int toaa = adjseg->aa_;
//			int tobb = adjseg->bb_;
//			int le = tobb - toaa;
//
//			if (seg->tube_id_ != adjseg->tube_id_)
//			{
//				// edge->seg->s_ = root->s_; s_ is not used
//				for (size_t i = 0; i < seg->colli_segs_.size(); i++)
//				{
//					if (seg->colli_segs_[i] == adjseg)
//					{
//						double c = seg->colli_ab_[i].x;
//						double d = seg->colli_ab_[i].y;
//
//						double rc = 0, rsc = 0, rd = 0, rsd = 0;
//
//						if (seg->b_ > seg->a_)
//						{
//							rc = double(c - seg->a_) / double(seg->b_ - seg->a_);
//							rd = double(d - seg->a_) / double(seg->b_ - seg->a_);
//						}
//
//						double cc = seg->aa_ + rc * (seg->bb_ - seg->aa_);
//						double dd = seg->aa_ + rd * (seg->bb_ - seg->aa_);
//
//						if (c == d)
//						{
//							adjseg->aa_ = int(cc - c + adjseg->a_);
//							adjseg->bb_ = int(cc - c + adjseg->b_);
//						}
//						else
//						{
//							double fac1 = dd - cc;
//							double fac2 = d - adjseg->a_;
//							double fac3 = d - c;
//							double fac4 = fac1 * fac2 / fac3;
//							adjseg->aa_ = int(dd - fac4);
//
//
//							fac2 = adjseg->b_ - c;
//							fac4 = fac1 * fac2 / fac3;
//							adjseg->bb_ = int(fac4 + cc);
//						}
//
//						if (adjseg->bb_ < adjseg->aa_)
//						{
//							std::cout << "Something strange in dfs_change_aabb" << std::endl;
//							adjseg->bb_ = adjseg->aa_;
//						}
//
//						break;
//					}
//				}
//			}
//			else
//			{
//				if (edge->direction == 0) // forward
//				{
//					if (adjseg->aa_ < oaa || adjseg->aa_ > obb + 1 ||
//						adjseg->bb_ < obb)
//					{
//						std::cout << "---------------------Wrong 1!" << std::endl;
//					}
//
//					if (adjseg->aa_ == obb + 1)	// special case
//					{
//						adjseg->aa_ = seg->bb_ + 1;
//						adjseg->bb_ = adjseg->aa_ + le;
//					}
//					else
//					{
//						if (oaa == obb)
//						{
//							if (oaa != adjseg->aa_)
//							{
//								std::cout << "ssss wrong!" << std::endl;
//							}
//							adjseg->aa_ = int((seg->aa_ + seg->bb_) / 2.0);
//							adjseg->bb_ = adjseg->aa_ + le;
//						}
//						else
//						{
//							double ratio = double(adjseg->aa_ - oaa) / double(obb - oaa);
//							adjseg->aa_ = int(seg->aa_ + (seg->bb_ - seg->aa_) * ratio);
//							adjseg->bb_ = adjseg->aa_ + le;
//						}
//					}
//				}
//				else
//				{
//					if (adjseg->bb_ < oaa - 1 || adjseg->bb_ > obb ||
//						adjseg->aa_ > oaa)
//					{
//						std::cout << "=====================Wrong 2!" << std::endl;
//					}
//
//					if (adjseg->bb_ == oaa - 1)
//					{
//						adjseg->bb_ = seg->aa_ - 1;
//						adjseg->aa_ = adjseg->bb_ - le;
//					}
//					else
//					{
//						if (oaa == obb)
//						{
//							if (adjseg->bb_ != obb)
//							{
//								std::cout << "SSSS wrong" << std::endl;
//							}
//							adjseg->bb_ = int((seg->aa_ + seg->bb_) / 2.0);
//							adjseg->aa_ = adjseg->bb_ - le;
//						}
//						else
//						{
//							double ratio = double(adjseg->bb_ - oaa) / double(obb - oaa);
//							adjseg->bb_ = int(seg->aa_ + (seg->bb_ - seg->aa_) * ratio);
//							adjseg->aa_ = adjseg->bb_ - le;
//						}
//					}
//				}
//
//			}
//			dfs_change_aabb(adjseg, toaa, tobb);
//		}
//		edge = edge->next;
//	}
//}

void change_pre_aabb(Segment * seg, Segment * preseg, int oaa, int obb, bool show)
{
	preseg->is_traversed = true;

	int le = preseg->bb_ - preseg->aa_;

	if (preseg->bb_ < oaa - 1 || preseg->bb_ > obb || preseg->aa_ > oaa)
	{
		std::cout << "=====================Wrong 2!" << std::endl;
		std::cout << "oaa:" << oaa << " obb:" << obb << " preseg->aa_:" << preseg->aa_ << " preseg->bb_:" << preseg->bb_ << std::endl;
		std::cout << "tubeid:" << seg->tube_id_ << std::endl;
	}

	if (preseg->bb_ == oaa - 1)
	{
		/*if (preseg->tube_id_ == 102)
		{
			std::cout << "change pre 102*******" << preseg->aa_ << " " << preseg->bb_ << " " << oaa << " " << seg->aa_ << std::endl;
		}*/
		preseg->bb_ = seg->aa_ - 1;
		preseg->aa_ = preseg->bb_ - le;
	}
	else
	{
		if (oaa == obb)
		{
			if (preseg->bb_ != obb)
			{
				std::cout << "SuBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB2 wrong" << std::endl;
			}

			preseg->bb_ = seg->aa_;
			preseg->aa_ = preseg->bb_ - le;
		}
		else
		{
			
			int rh = preseg->bb_ - oaa;
			preseg->bb_ = seg->aa_ + rh;
			if (preseg->bb_ > seg->bb_)
				preseg->bb_ = seg->bb_;
			preseg->aa_ = preseg->bb_ - le;
		}
	}
	
	if(show)
		std::cout << "pre:"<< preseg->aa_ << "," << preseg->bb_ << std::endl;

	/*if (preseg->aa_ < -100000 || preseg->aa_ > 100000 || preseg->bb_ < -100000 || preseg->bb_ > 100000)
	{

		std::cout << preseg->aa_ << "," << preseg->bb_ << "," << preseg->tube_id_ << "uuuuuuuuuuuuuuuuuu" << std::endl;
		std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
	}*/
}

void change_post_aabb(Segment* seg, Segment* postseg, int oaa, int obb, bool show)
{
	postseg->is_traversed = true;
	int le = postseg->bb_ - postseg->aa_;
	if (postseg->aa_ < oaa || postseg->aa_ > obb + 1 ||
		postseg->bb_ < obb)
	{
		std::cout << "---------------------Wrong 1!" << std::endl;
		std::cout << "oaa:" << oaa << " obb:" << obb << " postseg->aa_:" << postseg->aa_ << " postseg->bb_:" << postseg->bb_ << std::endl;
		std::cout << "tubeid:" << seg->tube_id_ << std::endl;
	}

	if (postseg->aa_ == obb + 1)	// special case
	{
		/*if (postseg->tube_id_ == 102)
		{
			std::cout << "change post 102*******" << postseg->aa_ << " " << postseg->bb_ << " " << oaa << " " << seg->aa_ << std::endl;
		}*/

		postseg->aa_ = seg->bb_ + 1;
		postseg->bb_ = postseg->aa_ + le;
	}
	else
	{
		if (oaa == obb)
		{
			if (oaa != postseg->aa_)
			{
				std::cout << "SuBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB1 wrong!" << std::endl;
			}
			postseg->aa_ = seg->bb_;
			postseg->bb_ = postseg->aa_ + le;
		}
		else
		{
			int lh = obb - postseg->aa_;
			postseg->aa_ = seg->bb_ - lh;
			if (postseg->aa_ < seg->aa_)
				postseg->aa_ = seg->aa_;
			postseg->bb_ = postseg->aa_ + le;
		}
	}

	if(show)
		std::cout << "post:"<< postseg->aa_ << "," << postseg->bb_ << std::endl;

	/*if (postseg->aa_ < -100000 || postseg->aa_ > 100000 || postseg->bb_ < -100000 || postseg->bb_ > 100000)
	{
		std::cout << postseg->aa_ << "," << postseg->bb_ << "," << postseg->tube_id_ << "vvvvvvvvvvvvvvvvvv" << std::endl;
	}*/
}

void change_colli_aabb(Segment* seg, Segment* adjseg, bool show)
{
	for (size_t i = 0; i < seg->colli_segs_.size(); i++)
	{
		if (seg->colli_segs_[i] == adjseg)
		{
			double c = seg->colli_ab_[i].x;
			double d = seg->colli_ab_[i].y;

			double rc = 0, rsc = 0, rd = 0, rsd = 0;

			if (seg->b_ > seg->a_)
			{
				rc = double(c - seg->a_) / double(seg->b_ - seg->a_);
				rd = double(d - seg->a_) / double(seg->b_ - seg->a_);
			}

			double cc = seg->aa_ + rc * (seg->bb_ - seg->aa_);
			double dd = seg->aa_ + rd * (seg->bb_ - seg->aa_);

			if (c == d)
			{
				adjseg->aa_ = int(cc - c + adjseg->a_);
				adjseg->bb_ = int(cc - c + adjseg->b_);
			}
			else
			{
				
				double fac1 = dd - cc;
				double fac2 = d - adjseg->a_;
				double fac3 = d - c;
				double fac4 = fac1 * fac2 / fac3;
				adjseg->aa_ = int(dd - fac4);

				/*if (adjseg->tube_id_ == 117)
				{
					std::cout << "aaaaa:" << rc << "," << rd << "," << cc << "," << dd << std::endl;
					std::cout << "bbbbb:" << seg->aa_ << "," << seg->bb_ << std::endl;
					std::cout << "ccccc:" << fac1 << "," << fac2 << "," << fac3 << "," << fac4 << ","<<adjseg->aa_<<std::endl;
				}*/

				fac2 = adjseg->b_ - c;
				fac4 = fac1 * fac2 / fac3;
				adjseg->bb_ = int(fac4 + cc);
			}

			if (adjseg->bb_ < adjseg->aa_)
			{
				std::cout << "Something strange in dfs_change_aabb" << std::endl;
				adjseg->bb_ = adjseg->aa_;
			}

			/*if (adjseg->aa_ < -100000 || adjseg->aa_ > 100000 || adjseg->bb_ < -100000 || adjseg->bb_ > 100000)
			{
				std::cout << adjseg->aa_ << "," << adjseg->bb_ << "," << adjseg->tube_id_ << "====================" << std::endl;
			}*/

			if(show)
				std::cout << "adj:" << adjseg->aa_ << "," << adjseg->bb_ << std::endl;

			break;
		}
	}
}

void check_consistency(Context & context)
{
	for (size_t i = 0; i < context.heads.size(); i++)
	{
		Segment* head = context.heads[i];
		
		head = head->next_;

		while (head != NULL)
		{
			if (head->next_ != NULL)
			{
				/*if (head->aa_ > head->next_->aa_ ||
					head->bb_ < head->next_->aa_ - 1 || head->bb_ > head->next_->bb_)*/
				if (head->bb_ != head->next_->aa_ - 1)
				{
					std::cout << "check consistency wrong! ................................." << std::endl;
					std::cout << "(" << head->tube_id_ << ": " << head->a_ << "->" << head->b_ << "," << head->aa_ << "->" << head->bb_ << ")" << std::endl;
					std::cout << "(" << head->next_->tube_id_ << ": " << head->next_->a_ << "->" << head->next_->b_ << "," << head->next_->aa_ << "->" << head->next_->bb_ << ")" << std::endl;
				}
			}
			head = head->next_;
		}
	}
}



void find_all_colli_segs(Segment* seg, Segment* tmp_segs[], int a[], int b[], int& p)
{
	seg->is_traversed = true;

	for (size_t i = 0; i < seg->colli_segs_.size(); i++)
	{
		if (!seg->colli_segs_[i]->is_traversed) {
			//tmp_segs.push_back(seg->colli_segs_[i]);
			tmp_segs[p] = seg->colli_segs_[i];
			int tmp_oaa = seg->colli_segs_[i]->aa_;
			int tmp_obb = seg->colli_segs_[i]->bb_;
			//tmp_aabb.push_back(cv::Point2i(tmp_oaa, tmp_obb));
			a[p] = tmp_oaa;
			b[p] = tmp_obb;
			p++;
			
			seg->colli_segs_[i]->aa_ = seg->aa_;
			seg->colli_segs_[i]->bb_ = seg->bb_;

			seg->colli_segs_[i]->scale_ = seg->scale_;
			
			if (seg->colli_segs_[i]->colli_segs_.size() != 0) {
				find_all_colli_segs(seg->colli_segs_[i], tmp_segs, a, b, p);
			}
		}

	}
}

void find_all_colli_segs(Segment* seg, std::vector<Segment*>& tmp_segs, std::vector<int>& a, std::vector<int>& b)
{
	seg->is_traversed = true;

	for (size_t i = 0; i < seg->colli_segs_.size(); i++)
	{
		if (!seg->colli_segs_[i]->is_traversed) {
			seg->colli_segs_[i]->is_traversed = true;
			tmp_segs.push_back(seg->colli_segs_[i]);

			int tmp_oaa = seg->colli_segs_[i]->aa_;
			int tmp_obb = seg->colli_segs_[i]->bb_;
			a.push_back(tmp_oaa);
			b.push_back(tmp_obb);

			/*if (seg->colli_segs_[i]->tube_id_ == 102)
			{
				std::cout << "change colli 102*******" << seg->colli_segs_[i]->aa_ << " " << seg->colli_segs_[i]->bb_ <<  " " <<" "<<seg->tube_id_<<" "<< seg->aa_ << " " << seg->bb_ << std::endl;
			}*/


			seg->colli_segs_[i]->aa_ = seg->aa_;
			seg->colli_segs_[i]->bb_ = seg->bb_;

			seg->colli_segs_[i]->scale_ = seg->scale_;

			/*if (seg->colli_segs_[i]->colli_segs_.size() != 0) {
				find_all_colli_segs(seg->colli_segs_[i], tmp_segs, a, b);
			}*/
		}

	}
}

void dfs_change_aabb(Segment* seg, int oaa, int obb, bool print_self)
{
	if (print_self && !seg->is_traversed)
	{
		std::cout << "(" << seg->tube_id_ << ": " << seg->a_ << "->" << seg->b_ << "," << seg->aa_ << "->" << seg->bb_ << ")" << std::endl;
	}

	seg->is_traversed = true;
/*
	Segment* tmp_segs[500];
	int a[500];
	int b[500];
	int p = 0;

	find_all_colli_segs(seg, tmp_segs, a, b, p);

	for (int i = 0; i < p; i++)
	{
		dfs_change_aabb(tmp_segs[i], a[i], b[i]);
	}*/

	std::vector<Segment*> tmp_segs;
	std::vector<int> a;
	std::vector<int> b;
	find_all_colli_segs(seg, tmp_segs, a, b);
	
	/*if (seg->tube_id_ == 77)
	{
		for (int i = 0; i < tmp_segs.size(); i++)
		{
			std::cout << "77 colli ssge" << tmp_segs[i]->tube_id_ << " " << tmp_segs[i]->a_ << " " << tmp_segs[i]->b_ << " " << tmp_segs[i]->aa_ << " " << tmp_segs[i]->bb_ << " " << std::endl;
		}
	}*/

	for (int i = 0; i < tmp_segs.size(); i++)
	{
		dfs_change_aabb(tmp_segs[i], a[i], b[i], print_self);
	}

	/*for (size_t i = 0; i < seg->colli_segs_.size(); i++)
	{
		if (!seg->colli_segs_[i]->is_traversed)
		{
			int tmp_oaa = seg->colli_segs_[i]->aa_;
			int tmp_obb = seg->colli_segs_[i]->bb_;
			seg->colli_segs_[i]->aa_ = seg->aa_;
			seg->colli_segs_[i]->bb_ = seg->bb_;
			dfs_change_aabb(seg->colli_segs_[i], tmp_oaa, tmp_obb, print_self);
		}
	}*/

	/*std::vector<cv::Point2i> tmp_aabb;
	std::vector<Segment*> tmp_segs;
	for (size_t i = 0; i < seg->colli_segs_.size(); i++)
	{
		if (!seg->colli_segs_[i]->is_traversed)
		{
			int tmp_oaa = seg->colli_segs_[i]->aa_;
			int tmp_obb = seg->colli_segs_[i]->bb_;
			tmp_aabb.push_back(cv::Point2i(tmp_oaa, tmp_obb));

			seg->colli_segs_[i]->aa_ = seg->aa_;
			seg->colli_segs_[i]->bb_ = seg->bb_;

			if (print_self)
			{
				std::cout << "(" << seg->colli_segs_[i]->tube_id_ << ": " << seg->colli_segs_[i]->a_ << "->" << seg->colli_segs_[i]->b_ << "," << seg->colli_segs_[i]->aa_ << "->" << seg->colli_segs_[i]->bb_ << ")" << std::endl;
			}

			seg->colli_segs_[i]->is_traversed = true;
			tmp_segs.push_back(seg->colli_segs_[i]);
		}
	}

	for (size_t i = 0; i < tmp_segs.size(); i++)
	{
		dfs_change_aabb(tmp_segs[i], tmp_aabb[i].x, tmp_aabb[i].y, print_self);
	}*/

	if (seg->next_ != NULL && !seg->next_->is_traversed)
	{
		int tmp_oaa = seg->next_->aa_;
		int tmp_obb = seg->next_->bb_;
		change_post_aabb(seg, seg->next_, oaa, obb, false);
		dfs_change_aabb(seg->next_, tmp_oaa, tmp_obb, print_self);
	}

	if (seg->last_ != NULL && !seg->last_->is_traversed)
	{
		int tmp_oaa = seg->last_->aa_;
		int tmp_obb = seg->last_->bb_;
		if (seg->last_->last_ != NULL) // not head
		{
			change_pre_aabb(seg, seg->last_, oaa, obb, false);
		}
		dfs_change_aabb(seg->last_, tmp_oaa, tmp_obb, print_self);
	}
}

void bfs_change_aabb(Context& context, Segment* seg, int oaa, int obb, bool show)
{
	seg->is_traversed = true;

	int tid = seg->tube_id_;
	int sid = 0;

	for (size_t i = 0; i < context.all_segs[tid].size(); i++)
	{
		if (context.all_segs[tid][i] == seg)
		{
			sid = int(i);
			break;
		}
	}

	int uoaa = oaa;
	int uobb = obb;
	for (int j = sid - 1; j >= 0; j--)
	{
		context.all_segs[tid][j]->is_traversed = true;
		int toaa = context.all_segs[tid][j]->aa_;
		int tobb = context.all_segs[tid][j]->bb_;
		change_pre_aabb(context.all_segs[tid][j + 1], context.all_segs[tid][j], uoaa, uobb, show);
		uoaa = toaa;
		uobb = tobb;
	}

	uoaa = oaa;
	uobb = obb;
	for (int j = sid + 1; j < int(context.all_segs[tid].size()); j++)
	{
		context.all_segs[tid][j]->is_traversed = true;
		int toaa = context.all_segs[tid][j]->aa_;
		int tobb = context.all_segs[tid][j]->bb_;
		change_post_aabb(context.all_segs[tid][j - 1], context.all_segs[tid][j], uoaa, uobb, show);
		uoaa = toaa;
		uobb = tobb;
	}

	for (size_t j = 0; j < context.all_segs[tid].size(); j++)
	{
	
		EdgeNode* edge = context.all_segs[tid][j]->firstedge;

		while (edge != NULL)
		{
			if (edge->is_retained && !edge->is_break
				&& !edge->seg->is_traversed
				&& edge->seg->tube_id_ != context.all_segs[tid][j]->tube_id_)
			{
				uoaa = edge->seg->aa_;
				uobb = edge->seg->bb_;
				change_colli_aabb(context.all_segs[tid][j], edge->seg, show);
				bfs_change_aabb(context, edge->seg, uoaa, uobb, show);
			}
			edge = edge->next;
		}
	}
}



//
//void dfs_change_aabb(Segment* seg, int oaa, int obb)
//{
//	seg->is_traversed = true;
//
//	std::vector<Segment*> self;
//	std::vector<Segment*> other;
//
//	EdgeNode* edge = seg->firstedge;
//
//	while (edge != NULL)
//	{
//		if (edge->is_retained && !edge->seg->is_traversed && !edge->is_break)
//		{
//			if (edge->seg->tube_id_ == seg->tube_id_)
//			{
//				self.push_back(edge->seg);
//			}
//			else
//			{
//				other.push_back(edge->seg);
//			}
//		}
//		edge = edge->next;
//	}
//
//	for (size_t i = 0; i < self.size(); i++)
//	{
//		Segment* adjseg = self[i];
//		int toaa = adjseg->aa_;
//		int tobb = adjseg->bb_;
//		int le = tobb - toaa;
//
//		if (adjseg->a_ > seg->a_) // forward
//		{
//			if (adjseg->aa_ < oaa || adjseg->aa_ > obb + 1 ||
//				adjseg->bb_ < obb)
//			{
//				std::cout << "---------------------Wrong 1!" << std::endl;
//			}
//
//			if (adjseg->aa_ == obb + 1)	// special case
//			{
//				adjseg->aa_ = seg->bb_ + 1;
//				adjseg->bb_ = adjseg->aa_ + le;
//			}
//			else
//			{
//				if (oaa == obb)
//				{
//					if (oaa != adjseg->aa_)
//					{
//						std::cout << "ssss wrong!" << std::endl;
//					}
//					adjseg->aa_ = int((seg->aa_ + seg->bb_) / 2.0);
//					adjseg->bb_ = adjseg->aa_ + le;
//				}
//				else
//				{
//					double ratio = double(adjseg->aa_ - oaa) / double(obb - oaa);
//					adjseg->aa_ = int(seg->aa_ + (seg->bb_ - seg->aa_) * ratio);
//					adjseg->bb_ = adjseg->aa_ + le;
//				}
//			}
//		}
//		else
//		{
//			if (adjseg->bb_ < oaa - 1 || adjseg->bb_ > obb ||
//				adjseg->aa_ > oaa)
//			{
//				std::cout << "=====================Wrong 2!" << std::endl;
//			}
//
//			if (adjseg->bb_ == oaa - 1)
//			{
//				adjseg->bb_ = seg->aa_ - 1;
//				adjseg->aa_ = adjseg->bb_ - le;
//			}
//			else
//			{
//				if (oaa == obb)
//				{
//					if (adjseg->bb_ != obb)
//					{
//						std::cout << "SSSS wrong" << std::endl;
//					}
//					adjseg->bb_ = int((seg->aa_ + seg->bb_) / 2.0);
//					adjseg->aa_ = adjseg->bb_ - le;
//				}
//				else
//				{
//					double ratio = double(adjseg->bb_ - oaa) / double(obb - oaa);
//					adjseg->bb_ = int(seg->aa_ + (seg->bb_ - seg->aa_) * ratio);
//					adjseg->aa_ = adjseg->bb_ - le;
//				}
//			}
//		}
//		dfs_change_aabb(adjseg, toaa, tobb);
//	}
//
//	for (size_t i = 0; i < other.size(); i++)
//	{
//		Segment* adjseg = other[i];
//		int toaa = adjseg->aa_;
//		int tobb = adjseg->bb_;
//		int le = tobb - toaa;
//		
//		for (size_t i = 0; i < seg->colli_segs_.size(); i++)
//		{
//			if (seg->colli_segs_[i] == adjseg)
//			{
//				double c = seg->colli_ab_[i].x;
//				double d = seg->colli_ab_[i].y;
//
//				double rc = 0, rsc = 0, rd = 0, rsd = 0;
//
//				if (seg->b_ > seg->a_)
//				{
//					rc = double(c - seg->a_) / double(seg->b_ - seg->a_);
//					rd = double(d - seg->a_) / double(seg->b_ - seg->a_);
//				}
//
//				double cc = seg->aa_ + rc * (seg->bb_ - seg->aa_);
//				double dd = seg->aa_ + rd * (seg->bb_ - seg->aa_);
//
//				if (c == d)
//				{
//					adjseg->aa_ = int(cc - c + adjseg->a_);
//					adjseg->bb_ = int(cc - c + adjseg->b_);
//				}
//				else
//				{
//					double fac1 = dd - cc;
//					double fac2 = d - adjseg->a_;
//					double fac3 = d - c;
//					double fac4 = fac1 * fac2 / fac3;
//					adjseg->aa_ = int(dd - fac4);
//
//
//					fac2 = adjseg->b_ - c;
//					fac4 = fac1 * fac2 / fac3;
//					adjseg->bb_ = int(fac4 + cc);
//				}
//
//				if (adjseg->bb_ < adjseg->aa_)
//				{
//					std::cout << "Something strange in dfs_change_aabb" << std::endl;
//					adjseg->bb_ = adjseg->aa_;
//				}
//
//				break;
//			}
//			dfs_change_aabb(adjseg, toaa, tobb);
//		}
//	}
//}

//
//void dfs_change_aabb(Segment* root)	// root's aabb has been changed
//{
//	root->is_traversed = true;
//	EdgeNode* edge = root->firstedge;
//	while (edge != NULL)	// transfer root's aabb to its neighbors
//	{
//		if (edge->is_retained && !edge->seg->is_traversed && !edge->is_break)
//		{
//			if (root->tube_id_ != edge->seg->tube_id_)
//			{
//				edge->seg->s_ = root->s_;
//				for (size_t i = 0; i < root->colli_segs_.size(); i++)
//				{
//					if (root->colli_segs_[i] == edge->seg)
//					{
//						/*
//						root       edge->seg
//
//									sa
//										|
//						a				| rsc
//						   |rc			|
//						c			c
//
//
//
//						d			d
//
//						b
//
//
//									sb
//
//						*/
//						double c = root->colli_ab_[i].x;
//						double d = root->colli_ab_[i].y;
//						double rc = 0, rsc = 0, rd = 0, rsd = 0;
//						if (root->b_ > root->a_)
//						{
//							rc = double(c - root->a_) / double(root->b_ - root->a_);
//							rd = double(d - root->a_) / double(root->b_ - root->a_);
//						}
//
//						if (edge->seg->b_ > edge->seg->a_)
//						{
//							rsc = double(c - edge->seg->a_) / double(edge->seg->b_ - edge->seg->a_);
//							rsd = double(d - edge->seg->a_) / double(edge->seg->b_ - edge->seg->a_);
//						}
//
//						double cc = root->aa_ + rc * (root->bb_ - root->aa_);
//						double dd = root->aa_ + rd * (root->bb_ - root->aa_);
//						double slen = edge->seg->len_ / edge->seg->s_;
//
//						edge->seg->aa_ = int(cc - rsc * (slen - 1) + 0.5);
//						edge->seg->bb_ = int(edge->seg->aa_ + slen - 0.5);
//						if (edge->seg->bb_ < edge->seg->aa_)
//						{
//							edge->seg->bb_ = edge->seg->aa_;
//						}
//
//						break;
//					}
//				}
//			}
//			else
//			{
//				if (edge->direction == 0)
//				{
//					edge->seg->aa_ = root->bb_ + 1;
//					edge->seg->bb_ = int(edge->seg->len_ / edge->seg->s_ + edge->seg->aa_ - 1);
//					if (edge->seg->bb_ < edge->seg->aa_)
//					{
//						edge->seg->bb_ = edge->seg->aa_;
//					}
//
//				}
//				else
//				{
//					edge->seg->bb_ = root->aa_ - 1;
//					edge->seg->aa_ = int(edge->seg->bb_ - edge->seg->len_ / edge->seg->s_ + 1);
//					if (edge->seg->bb_ < edge->seg->aa_)
//					{
//						edge->seg->bb_ = edge->seg->aa_;
//					}
//				}
//
//			}
//			dfs_change_aabb(edge->seg);
//		}
//		edge = edge->next;
//	}
//}


cv::Mat create_segment_graph(std::string filename, std::vector<std::vector<Segment*>>& all_segs,
	std::vector<std::pair<EdgeNode*, EdgeNode*>>& all_edges,
	int seg_len,
	int & srclen)
{
	std::ifstream infile(filename);
	std::string line;

	std::getline(infile, line);
	int number1;
	int max_frame;
	std::istringstream(line) >> number1;				// number of tubes

	cv::Mat ret = cv::Mat(number1, number1, CV_8UC1, cv::Scalar(0));
	
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
		std::istringstream(line) >> a >> b >> c >> d >> e >> f;	// a: tube1, d: tube2, (b,c): the collision segment
		tube_colli_pre[a - 1].push_back(cv::Point3i(b - 1, c, d));
		ret.at<uchar>(a - 1, b - 1) = 1;
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
				if (tube_colli[i][k].x == id)     // combine two segments of the same tube
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


#ifdef __DEBUG__

	std::cout << "Should be no multiple correspondence between two tubes" << std::endl;
	for (size_t i = 0; i < tube_colli.size(); i++)
	{
		for (size_t j = 0; j < tube_colli[i].size(); j++)
		{
			std::cout << i << "," << tube_colli[i][j].x << "," << tube_colli[i][j].y << "," << tube_colli[i][j].z << std::endl;
		}
	}

#endif // __DEBUG__

#ifdef __DEBUG__
	for (size_t i = 0; i < tube_colli.size(); i++)
	{
		std::cout << "tube:" << i << std::endl;
		for (size_t j = 0; j < tube_colli[i].size(); j++)
		{
			std::cout << tube_colli[i][j] << std::endl;
		}
	}

	std::cout << "##########################################" << std::endl;
#endif

	std::vector<std::vector<cv::Point2i>> combined;			// a tube may collide with many other tubes, therefore the collision points with different tubes should be merged
	combined.resize(number1);
	std::vector<std::vector<std::vector<int>>> combined_idx;	// combined_idx remember the tubes that collide with current tube
	combined_idx.resize(number1);

	for (int i = 0; i < number1; i++)
	{
		std::vector<std::vector<int>> flag;	// record collided tubes at each frame
		flag.resize(max_frame);

		for (size_t j = 0; j < tube_colli[i].size(); j++)
		{
			int n, a, b;
			n = tube_colli[i][j].x;
			a = tube_colli[i][j].y;
			b = tube_colli[i][j].z;

			for (int k = a; k <= b; k++)
			{
				flag[k].push_back(n);
			}
		}

		int beg = -1;
		for (size_t j = 0; j < flag.size(); j++)
		{
			if (!flag[j].empty() && beg == -1)
			{
				beg = int(j);
			}

			if (flag[j].empty() && beg != -1)
			{
				combined[i].push_back(cv::Point2i(beg, int(j) - 1));
				beg = -1;
			}
		}

		if (beg != -1)
		{
			combined[i].push_back(cv::Point2i(beg, int(flag.size()) - 1));
		}

		for (size_t j = 0; j < combined[i].size(); j++)
		{
			std::set<int> tmp;
			std::vector<int> tmp2;

			int a = combined[i][j].x;
			int b = combined[i][j].y;

			for (int k = a; k <= b; k++)
			{
				for (size_t kk = 0; kk < flag[k].size(); kk++)
				{
					tmp.insert(flag[k][kk]);
				}
			}

			std::set<int>::iterator iter;
			for (iter = tmp.begin(); iter != tmp.end(); iter++)
			{
				tmp2.push_back(*iter);
			}
			combined_idx[i].push_back(tmp2);
		}
	}

#ifdef __DEBUG__
	for (size_t i = 0; i < combined.size(); i++)
	{
		std::cout << "tube:" << i << std::endl;

		for (size_t j = 0; j < combined[i].size(); j++)
		{
			std::cout << combined[i][j] << std::endl;
			for (size_t k = 0; k < combined_idx[i][j].size(); k++)
			{
				std::cout << combined_idx[i][j][k] << ",";
			}
			std::cout << std::endl;
		}
	}
#endif

	for (size_t i = 0; i < tube_ab.size(); i++)					// Get all segments of all tubes
	{
		int A = tube_ab[i].x;
		int B = tube_ab[i].y;

		std::vector<Segment*> segs;

		if (combined[i].empty())
		{
			extract_segment_from(A, B, int(i), seg_len, segs);
		}
		else
		{
			// before all collision segments
			if (A < combined[i][0].x)
			{
				extract_segment_from(A, combined[i][0].x - 1, int(i), seg_len, segs);
			}
			Segment* seg = new Segment(combined[i][0].x, combined[i][0].y, int(i),
				combined[i][0].y - combined[i][0].x + 1, true);

			segs.push_back(seg);

			for (int j = 1; j < int(combined[i].size()); j++)
			{
				extract_segment_from(combined[i][j - 1].y + 1, combined[i][j].x - 1,
					int(i), seg_len, segs);

				Segment* seg = new Segment(combined[i][j].x, combined[i][j].y, int(i),
					combined[i][j].y - combined[i][j].x + 1, true);

				segs.push_back(seg);
			}

			if (combined[i][combined[i].size() - 1].y + 1 <= B)
			{
				extract_segment_from(combined[i][combined[i].size() - 1].y + 1, B, int(i), seg_len, segs);
			}
		}

		all_segs.push_back(segs);
	}

#ifdef __DEBUG__
	for (size_t i = 0; i < all_segs.size(); i++)
	{
		for (size_t j = 0; j < all_segs[i].size(); j++)
		{
			std::cout << *(all_segs[i][j]) << std::endl;
		}
	}
#endif

	for (size_t i = 0; i < int(all_segs.size()); i++)
	{
		for (int j = 0; j<int(all_segs[i].size()) - 1; j++)
		{
			EdgeNode* edge1 = new EdgeNode;
			EdgeNode* edge2 = new EdgeNode;

			edge1->pre_seg = all_segs[i][j];
			edge1->seg = all_segs[i][j + 1];
			edge1->next = NULL;
			edge1->is_retained = false;
			edge1->is_traversed = false;
			edge1->is_break = false;
			edge1->direction = 0;
			edge1->weight = (all_segs[i][j]->len_ + all_segs[i][j + 1]->len_) / 2.0;

			edge2->pre_seg = all_segs[i][j + 1];
			edge2->seg = all_segs[i][j];
			edge2->next = NULL;
			edge2->is_retained = false;
			edge2->is_traversed = false;
			edge2->is_break = false;
			edge2->direction = 1;
			edge2->weight = (all_segs[i][j]->len_ + all_segs[i][j + 1]->len_) / 2.0;

			if (all_segs[i][j]->firstedge == NULL)
			{
				all_segs[i][j]->firstedge = edge1;
			}
			else
			{
				EdgeNode* node = all_segs[i][j]->firstedge;
				while (node->next != NULL)
				{
					node = node->next;
				}
				node->next = edge1;
			}

			if (all_segs[i][j + 1]->firstedge == NULL)
			{
				all_segs[i][j + 1]->firstedge = edge2;
			}
			else
			{
				EdgeNode* node = all_segs[i][j + 1]->firstedge;
				while (node->next != NULL)
				{
					node = node->next;
					//node->next = edge2;  // bug? but this else can't be reached
				}
				node->next = edge2;
			}

			all_edges.push_back(std::pair<EdgeNode*, EdgeNode*>(edge1, edge2));
		}
	}

	for (size_t i = 0; i < tube_colli.size(); i++)
	{
		for (size_t j = 0; j < tube_colli[i].size(); j++)
		{
			int id = tube_colli[i][j].x;
			int aa = tube_colli[i][j].y;
			int bb = tube_colli[i][j].z;

			Segment* seg1 = NULL, * seg2 = NULL;

			for (size_t k = 0; k < all_segs[i].size(); k++)
			{
				if (all_segs[i][k]->colli_ &&
					aa >= all_segs[i][k]->a_ &&
					bb <= all_segs[i][k]->b_)
				{
					seg1 = all_segs[i][k];
					break;
				}
			}

			for (size_t k = 0; k < all_segs[id].size(); k++)
			{
				if (all_segs[id][k]->colli_ &&
					aa >= all_segs[id][k]->a_ &&
					bb <= all_segs[id][k]->b_)
				{
					seg2 = all_segs[id][k];
					break;
				}
			}

			EdgeNode* en1 = seg1->firstedge;
			EdgeNode* en2 = seg2->firstedge;

			bool ff = false;

			while (en1 != NULL && en1->next != NULL)
			{
				en1 = en1->next;
				if (en1->seg == seg2)
				{
					ff = true;
				}
			}

			while (en2 != NULL && en2->next != NULL)
			{
				en2 = en2->next;
				if (en2->seg == seg1)
				{
					ff = true;
				}
			}

			if (!ff)
			{
				EdgeNode* edge1 = new EdgeNode;
				EdgeNode* edge2 = new EdgeNode;


				edge1->pre_seg = seg1;
				edge1->seg = seg2;
				edge1->weight = 1e10;
				edge1->next = NULL;
				edge1->is_retained = false;
				edge1->is_traversed = false;
				edge1->is_break = false;
				if (en1 == NULL)
				{
					seg1->firstedge = edge1;
				}
				else
				{
					en1->next = edge1;
				}

				edge2->pre_seg = seg2;
				edge2->seg = seg1;
				edge2->weight = 1e10;
				edge2->next = NULL;
				edge2->is_retained = false;
				edge2->is_traversed = false;
				edge2->is_break = false;
				if (en2 == NULL)
				{
					seg2->firstedge = edge2;
				}
				else
				{
					en2->next = edge2;
				}

				seg1->colli_segs_.push_back(seg2);
				seg1->colli_ab_.push_back(cv::Point2i(aa, bb));

				seg2->colli_segs_.push_back(seg1);
				seg2->colli_ab_.push_back(cv::Point2i(aa, bb));

				all_edges.push_back(std::pair<EdgeNode*, EdgeNode*>(edge1, edge2));
			}
		}
	}

	// Kruskal Algorithm

#ifdef __DEBUG__
	std::cout << "Kruskal:" << std::endl;
#endif

	for (size_t i = 0; i < all_edges.size(); i++)
	{
		int max_idx = -1;
		double max_value = -1e10;

		for (size_t j = 0; j < all_edges.size(); j++)
		{
			if (all_edges[j].first->is_traversed == false &&
				all_edges[j].first->weight > max_value)
			{
				max_value = all_edges[j].first->weight;
				max_idx = (int)j;
			}
		}

		all_edges[max_idx].first->is_traversed = true;
		all_edges[max_idx].second->is_traversed = true;

		reset_traverse(all_segs);

		if (check_connection(all_edges[max_idx].first->pre_seg, all_edges[max_idx].first->seg))
		{

		}
		else
		{
			all_edges[max_idx].first->is_retained = true;	// retain is important, traverse is just for programming
			all_edges[max_idx].second->is_retained = true;

#ifdef __DEBUG__
			std::cout << *(all_edges[max_idx].first->pre_seg);
			std::cout << *(all_edges[max_idx].first->seg);
			std::cout << std::endl;
#endif 
		}
	}

#ifdef __DEBUG__
	std::cout << "Print Retained Edges:" << std::endl;
	for (size_t i = 0; i < all_edges.size(); i++)
	{
		EdgeNode* edge1 = all_edges[i].first;
		EdgeNode* edge2 = all_edges[i].second;

		if (edge1->is_retained && edge2->is_retained)
		{
			std::cout << *(edge1->pre_seg);
			std::cout << *(edge1->seg);
			std::cout << edge1->weight << std::endl;
			std::cout << std::endl;
		}
	}
#endif
	return ret;
}

void test_extract_segment_from()
{
	std::vector<Segment*> segs;
	extract_segment_from(4, 3, 1, 100, segs);
	for (size_t i = 0; i < segs.size(); i++)
	{
		std::cout << *(segs[i]);
	}

	for (size_t i = 0; i < segs.size(); i++)
	{
		delete segs[i];
	}
}

