/*
 Copyright (C) Intel Corp.  2006.  All Rights Reserved.
 Intel funded Tungsten Graphics (http://www.tungstengraphics.com) to
 develop this 3D driver.

 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:

 The above copyright notice and this permission notice (including the
 next paragraph) shall be included in all copies or substantial
 portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE COPYRIGHT OWNER(S) AND/OR ITS SUPPLIERS BE
 LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 **********************************************************************/
 /*
  * Authors:
  *   Keith Whitwell <keith@tungstengraphics.com>
  */


#ifndef BRW_VS_H
#define BRW_VS_H


#include "brw_context.h"
#include "brw_eu.h"


struct brw_vs_prog_key {
   unsigned program_string_id;
   unsigned nr_userclip:4;
   unsigned copy_edgeflag:1;
   unsigned know_w_is_one:1;
   unsigned pad:26;
};


struct brw_vs_compile {
   struct brw_compile func;
   struct brw_vs_prog_key key;
   struct brw_vs_prog_data prog_data;

   const struct brw_vertex_program *vp;

   unsigned nr_inputs;

   unsigned first_output;
   unsigned nr_outputs;

   unsigned first_tmp;
   unsigned last_tmp;

   struct brw_reg r0;
   struct brw_reg r1;
   struct brw_reg regs[12][128];
   struct brw_reg tmp;
   struct brw_reg stack;

   struct {
       boolean used_in_src;
       struct brw_reg reg;
   } output_regs[128];

   struct brw_reg userplane[6];

};

void brw_vs_emit( struct brw_vs_compile *c );

#endif
