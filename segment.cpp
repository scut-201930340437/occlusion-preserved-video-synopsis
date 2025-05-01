#include "segment.h"
#include <iostream>  
using namespace std;

void split_segment_using_a(Segment* seg, int pos, Segment*& seg1, Segment*& seg2)
{
	int a = seg->a_;
	int b = seg->b_;
	int id = seg->tube_id_;

	seg1 = new Segment(a, pos - 1, id, pos - a);
	seg2 = new Segment(pos, b, id, b - pos + 1);
	
	seg1->circulation_id_ = seg->circulation_id_;
	seg2->circulation_id_ = seg->circulation_id_;
}

void split_segment_using_b(Segment* seg, int pos, Segment*& seg1, Segment*& seg2)
{
	int a = seg->a_;
	int b = seg->b_;
	int id = seg->tube_id_;

	seg1 = new Segment(a, pos, id, pos - a + 1);
	seg2 = new Segment(pos + 1, b, id, b - pos);

	seg1->circulation_id_ = seg->circulation_id_;
	seg2->circulation_id_ = seg->circulation_id_;
}

Segment::Segment(int a, int b, int n, int len, bool colli)
{
	a_ = a;
	b_ = b;
	aa_ = a;
	bb_ = b;
	tube_id_ = n;
	len_ = len;
	//s_ = 1;
	scale_ = 1.0;
	firstedge = NULL;
	colli_ = colli;
	is_traversed = false;
	circulation_traversed = false;
	//circulation_id_ = -1;
	//flowin_seg_ = NULL;
	next_ = NULL;
	last_ = NULL;
	is_next_break_ = false;
	is_last_break_ = false;

	divide_seg1_ = NULL;
	divide_seg2_ = NULL;

	seg_index_ = -1;
}


Segment::Segment(int a, int b, int newBeg, int newEnd, int n, int len)
{
	a_ = a;
	b_ = b;
	aa_ = newBeg;
	bb_ = newEnd;
	tube_id_ = n;
	len_ = len;
	//s_ = 1;
	scale_ = 1.0;
	firstedge = NULL;
	colli_ = false;
	is_traversed = false;
	circulation_traversed = false;
	//circulation_id_ = -1;
	//flowin_seg_ = NULL;
	next_ = NULL;
	last_ = NULL;
	is_next_break_ = false;
	is_last_break_ = false;

	divide_seg1_ = NULL;
	divide_seg2_ = NULL;

	seg_index_ = -1;
}


Segment::Segment(const Segment& seg)
{
	this->a_ = seg.a_;
	this->b_ = seg.b_;
	this->aa_ = seg.aa_;
	this->bb_ = seg.bb_;
	//this->s_ = seg.s_;
	this->len_ = seg.len_;
	this->tube_id_ = seg.tube_id_;
	this->colli_ = seg.colli_;
	this->is_traversed = seg.is_traversed;
	this->circulation_traversed = seg.circulation_traversed;
	this->circulation_id_= seg.circulation_id_;
	this->firstedge = NULL;
	this->next_ = NULL;
	this->last_ = NULL;
	this->is_next_break_ = seg.is_next_break_;
	this->is_last_break_ = seg.is_last_break_;

	this->divide_seg1_ = seg.divide_seg1_;
	this->divide_seg2_ = seg.divide_seg2_;

	this->seg_index_ = seg.seg_index_;
	this->scale_ = seg.scale_;
	//this->flowin_seg_ = NULL;

	/// @todo: other members to be added
}


Segment::~Segment()
{

}

std::ostream& operator<<(std::ostream& out, const Segment& seg)
{
	out << "Begin:" << seg.a_ << ","
		<< "End:" << seg.b_ << ","
		<< "Length:" << seg.len_ << ","
		<< "Tube:" << seg.tube_id_ << '\n';
	return out;
}
//
//void Segment::change_accordingly(Segment* last_seg, int loaa, int lobb)
//{
//	if (last_seg != this->flowin_seg_)
//	{
//		std::cout << "+++++++++++++++++++++Wrong 3!" << std::endl;
//	}
//	int oaa = this->aa_;
//	int obb = this->bb_;
//
//	if (last_seg->tube_id_ != this->tube_id_)
//	{
//		this->s_ = last_seg->s_;
//		for (size_t i = 0; i < last_seg->colli_segs_.size(); i++)
//		{
//			if (last_seg->colli_segs_[i] == this)
//			{
//				/*
//				root       edge->seg
//
//							sa
//								|
//				a				| rsc
//				   |rc			|
//				c			c
//
//
//
//				d			d
//
//				b
//
//
//							sb
//
//				*/
//				double c = last_seg->colli_ab_[i].x;
//				double d = last_seg->colli_ab_[i].y;
//				double rc = 0, rsc = 0, rd = 0, rsd = 0;
//				
//				if (last_seg->b_ > last_seg->a_)
//				{
//					rc = double(c - last_seg->a_) / double(last_seg->b_ - last_seg->a_);
//					rd = double(d - last_seg->a_) / double(last_seg->b_ - last_seg->a_);
//				}
//
//				if (this->b_ > this->a_)
//				{
//					rsc = double(c - this->a_) / double(this->b_ - this->a_);
//					rsd = double(d - this->a_) / double(this->b_ - this->a_);
//				}
//
//				double cc = last_seg->aa_ + rc * (last_seg->bb_ - last_seg->aa_);
//				double dd = last_seg->aa_ + rd * (last_seg->bb_ - last_seg->aa_);
//				double slen = this->len_ / this->s_;
//				
//				this->aa_ = int(cc - rsc * (slen - 1) + 0.5);
//				this->bb_ = int(this->aa_ + slen - 0.5);
//				if (this->bb_ < this->aa_)
//				{
//					this->bb_ = this->aa_;
//				}
//
//				break;
//			}
//		}
//	}
//	else
//	{
//		if (direction() == 1) // last is preceding 
//		{
//			if (this->aa_ < loaa || this->aa_ > lobb + 1)
//			{
//				std::cout << "---------------------Wrong 1!" << std::endl;
//			}
//			
//			if (this->aa_ == lobb + 1)	// special case
//			{
//				int le = this->bb_ - this->aa_;
//				this->aa_ = last_seg->bb_ + 1;
//				this->bb_ = this->aa_ + le;
//			}
//			else
//			{
//				int le = this->bb_ - this->aa_;
//				double ratio = double(this->aa_ - loaa) / double(lobb - loaa);
//				this->aa_ = int(last_seg->aa_ + (last_seg->bb_ - last_seg->aa_) * ratio);
//				this->bb_ = this->aa_ + le;
//			}
//		}
//		else if(direction() == 0) // last is from back
//		{
//			if (this->bb_ < loaa - 1 || this->bb_ > lobb)
//			{
//				std::cout << "=====================Wrong 2!" << std::endl;
//			}
//
//			if (this->bb_ = loaa - 1)
//			{
//				int le = this->bb_ - this->aa_;
//				this->bb_ = last_seg->aa_ - 1;
//				this->aa_ = this->bb_ - le;
//			}
//			else
//			{
//				int le = this->bb_ - this->aa_;
//				double ratio = double(this->bb_ - loaa) / double(lobb - loaa);
//				this->bb_ = int(last_seg->aa_ + (last_seg->bb_ - last_seg->aa_) * ratio);
//				this->aa_ = this->bb_ - le;
//			}
//			
//		}
//	}
//
//	for (size_t i = 0; i < flowout_segs_.size(); i++)
//	{
//		flowout_segs_[i]->change_accordingly(this, oaa, obb);
//	}
//}

void Segment::split_using_a(int pos, Segment * s1, Segment * s2)
{
	Segment* seg1, * seg2;
	split_segment_using_a(this, pos, seg1, seg2);

	seg1->next_ = seg2;
	seg2->last_ = seg1;

	seg1->last_ = this->last_;
	seg2->next_ = this->next_;

	this->last_->next_ = seg1;
	if (this->next_ != NULL)
	{
		this->next_->last_ = seg2;
	}

	seg1->colli_segs_.push_back(s1);
	seg2->colli_segs_.push_back(s2);

	s1->colli_segs_.push_back(seg1);
	s2->colli_segs_.push_back(seg2);

	this->is_traversed = true;

	for (size_t i = 0; i < this->colli_segs_.size(); i++)
	{
		if (!this->colli_segs_[i]->is_traversed)
		{
			this->colli_segs_[i]->split_using_a(pos, seg1, seg2);
			delete this->colli_segs_[i];
		}
	}
}

void Segment::split_using_b(int pos, Segment* s1, Segment* s2)
{
	Segment* seg1, * seg2;
	split_segment_using_b(this, pos, seg1, seg2);

	seg1->next_ = seg2;
	seg2->last_ = seg1;

	seg1->last_ = this->last_;
	seg2->next_ = this->next_;

	this->last_->next_ = seg1;
	if (this->next_ != NULL)
	{
		this->next_->last_ = seg2;
	}

	seg1->colli_segs_.push_back(s1);
	seg2->colli_segs_.push_back(s2);

	s1->colli_segs_.push_back(seg1);
	s2->colli_segs_.push_back(seg2);

	this->is_traversed = true;

	for (size_t i = 0; i < this->colli_segs_.size(); i++)
	{
		if (!this->colli_segs_[i]->is_traversed)
		{
			this->colli_segs_[i]->split_using_b(pos, seg1, seg2);
			delete this->colli_segs_[i];
		}
	}
}

void Segment::collect_right_neighbors(std::vector<Segment*>& ans)
{
	this->is_traversed = true;
	if (this->next_ != NULL)
	{
		ans.push_back(this->next_);
	}

	for (size_t i = 0; i < colli_segs_.size(); i++)
	{
		if (!colli_segs_[i]->is_traversed)
		{
			colli_segs_[i]->collect_right_neighbors(ans);
		}
	}
}

void Segment::collect_left_neighbors(std::vector<Segment*>& ans)
{
	this->is_traversed = true;
	if (this->last_ != NULL)
	{
		ans.push_back(this->last_);
	}

	for (size_t i = 0; i < colli_segs_.size(); i++)
	{
		if (!colli_segs_[i]->is_traversed)
		{
			colli_segs_[i]->collect_left_neighbors(ans);
		}
	}
}