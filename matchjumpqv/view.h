/* view.h */

typedef struct view_{
  char* scaffold;
  int64_t view_start_addr;
  int64_t view_end_addr;
  int32_t view_start_scaffold_addr;
  int32_t view_end_scaffold_addr;
  uint32_t image_width;
  uint32_t image_hight;
  int unit;
  int flipimage;
  double xconversion_const;
  double yconversion_const;
}view;
