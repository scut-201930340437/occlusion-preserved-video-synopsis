#include "defines.h"
#include "segment.h"

void release_context(Context& context)
{
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			delete context.all_segs[i][j];
		}
		context.all_segs[i].clear();
	}
	context.all_segs.clear();

	for (size_t i = 0; i < context.all_edges.size(); i++)
	{
		delete context.all_edges[i].first;
		delete context.all_edges[i].second;
	}
	context.all_edges.clear();
}