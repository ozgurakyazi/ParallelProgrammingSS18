#ifndef _X_GRADIENT_PAR
#define _X_GRADIENT_PAR

#include <boost/gil/extension/io/jpeg_dynamic_io.hpp>
#include <omp.h>
#include"x_gradient.h"
#include <mutex>
#include <future>
using namespace boost::gil;
using namespace std;

template <typename Out> struct halfdiff_cast_channels; // forward declaration

const int step_size = 10; // 5 rows at a time
int curr_line = 0;
int total_rows;
int total_cols;
mutex curr_line_mutex;

int access_curr_line(){
  lock_guard<mutex> guard(curr_line_mutex);
  int val = curr_line;
  curr_line +=step_size;
  return val;
}

int get_new_row(){
  int my_row_start = access_curr_line();
  if(my_row_start>=total_rows){
    return -1;
  }
  return my_row_start;
}

template <typename SrcView, typename DstView>
void thread_work(const SrcView& src, const DstView& dst){
  int my_row_start = get_new_row();
  int my_row_end = min(my_row_start+step_size,total_rows);
  typedef typename channel_type<DstView>::type dst_channel_t;

  while(my_row_start != -1 ){
    while(my_row_start<my_row_end){
      typename SrcView::x_iterator src_it = src.row_begin(my_row_start);
      typename DstView::x_iterator dst_it = dst.row_begin(my_row_start);

      for (int x=1; x<total_cols; ++x) {
          static_transform(src_it[x-1], src_it[x+1], dst_it[x],halfdiff_cast_channels<dst_channel_t>());
      }
      my_row_start++;
    }
    my_row_start = get_new_row();
    my_row_end = min(my_row_start+step_size,total_rows);

  }
}


template <typename SrcView, typename DstView>
void x_gradient(const SrcView& src, const DstView& dst, int num_threads) {

  typedef typename channel_type<DstView>::type dst_channel_t;
  //cout<<src.height()<<endl;
  vector<future<void>> fut_vec;
  total_rows = src.height();
  total_cols = src.width() -1 ;


  for (size_t i = 0; i < num_threads; i++) {
      fut_vec.push_back(async(launch::async, [src, dst](){
          thread_work(src,dst);
      }));
  }

  for (auto &e: fut_vec){
    e.get();
  }

}





#endif // !_X_GRADIENT_PAR_
