/* uvid_decompress.cpp
   CSC 485B/578B/SENG 480B - Data Compression - Summer 2020

   Starter code for Assignment 5
   
   This placeholder code reads the (basically uncompressed) data produced by
   the uvid_compress starter code and outputs it in the uncompressed 
   YCbCr (YUV) format used for the sample video input files. To play the 
   the decompressed data stream directly, you can pipe the output of this
   program to the ffplay program, with a command like 

     ffplay -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 352x288 - 2>/dev/null
   (where the resolution is explicitly given as an argument to ffplay).

   B. Bird - 07/15/2020
*/

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include <tuple>
#include <list>
#include <math.h> 
#include "input_stream.hpp"
#include "yuv_stream.hpp"
#include "uv_common.hpp"


InputBitStream input_stream {std::cin};
int mblock_size = 16;



std::vector<std::vector<double>> Reverse_ZigZagBlock(std::vector<double> block_slice){

    std::vector<std::vector<double>> block = create_2d_vector<double>(8,8);
    std::list<double> block_slice_list(block_slice.begin(),block_slice.end());

    int row =0;
    int column =0;
    int inc = 1;

    //zig zag till half ( longest diagonal starting from index 7,0)
    for(unsigned int i =0; i<block.size(); i++){
        row =(i%2==0)?0:i;
        column = (i%2==0)?i:0;
        inc =(i%2==0)?-1:1;

        block.at(row).at(column) = block_slice_list.front();
        block_slice_list.pop_front();
        for(unsigned int j =0; j<i; j++){
            row+=-inc;
            column+=inc;
            block.at(row).at(column) = block_slice_list.front();
            block_slice_list.pop_front();
        }
    }

    //zig from longest dianoal to end ( index 7,7)
    for(unsigned int i =1; i<block.size(); i++){
        row =(i%2==0)?7:i;
        column = (i%2==0)?i:7;
        inc =(i%2==0)?1:-1;


        for(unsigned int j =i; j<block.size(); j++){
            block.at(row).at(column) = block_slice_list.front();
            block_slice_list.pop_front();
            row+=-inc;
            column+=inc;
        }
        
    }

    return block;


}

std::vector<std::vector<double>> Reverse_ZigZagOrder(std::vector<double> zig_zag, unsigned int height, unsigned int width){

    std::vector<std::vector<double>> plane = create_2d_vector<double>(height,width);
    auto block = create_2d_vector<double>(8,8);

    int slice_start = 0;
    for(unsigned int y = 0; y < height; y+=8){
        for (unsigned int x = 0; x < width; x+=8){
            
            
            auto begin = zig_zag.begin() + slice_start;
            auto end = zig_zag.begin() + slice_start + 64;

            std::vector<double> block_slice(begin, end); 
            
           

            auto ZigZag_block = Reverse_ZigZagBlock(block_slice);

            for(unsigned int i = 0; i < 8; i++){
                for(unsigned int j = 0; j < 8; j++){
                    if(((y+i)>=height || (x+j) >=width))
                        break;
                    plane.at(y+i).at(x+j)=ZigZag_block.at(i).at(j);
                }
            }
            slice_start+=64;
        }
    }

    return plane;

}


int read_variable_bits(){

    unsigned int sign = input_stream.read_bit();

    unsigned int bits = 0;
    while (input_stream.read_bit())
        bits++;

    int val = input_stream.read_bits(bits);
    val = (sign==0)?val:(val*-1);


    return val;
}

int main(int argc, char** argv){
 
    

    u32 height {input_stream.read_u32()};
    u32 width {input_stream.read_u32()};

    unsigned int quantization_factor = input_stream.read_bits(2);

    double q_factor = 0;
    switch(quantization_factor){
        case 1:
            q_factor = 2;
            break;
        case 2:
            q_factor = 1;
            break;
        case 3:
            q_factor = 0.5;
            break;
    }

    YUVStreamWriter writer {std::cout, width, height};

    auto last_Y = create_2d_vector<double>(height,width);
    auto last_Cb = create_2d_vector<double>(height/2,width/2);
    auto last_Cr = create_2d_vector<double>(height/2,width/2);
    
    unsigned int frame_num =0;

    while (input_stream.read_byte()){
        YUVFrame420& frame = writer.frame();
        frame_num++;

        unsigned int frame_type = input_stream.read_bits(2);

        std::vector<std::vector<int>> motion_vectors;
        std::vector<int> motion_v;
        //motion vectors
        if(frame_type){
            unsigned int motion_vector_size = read_variable_bits();
            
            for(int i =0; i<motion_vector_size; i++){
                for(int j =0; j<4; j++){
                    motion_v.push_back(read_variable_bits());
                }
                motion_vectors.push_back(motion_v);
                motion_v.clear();
            }
        }

        unsigned int ZigZag_Y_size = read_variable_bits();
        unsigned int ZigZag_Cb_size = read_variable_bits();
        unsigned int ZigZag_Cr_size = read_variable_bits();
        

        auto Y = create_2d_vector<double>(height,width);
        auto Cb_scaled = create_2d_vector<double>((height+1)/2,(width+1)/2);
        auto Cr_scaled = create_2d_vector<double>((height+1)/2,(width+1)/2);


        std::vector<double> ZigZag_Y;
        std::vector<double> ZigZag_Cb;
        std::vector<double> ZigZag_Cr;

        std::ofstream debugfile;
        debugfile.open ("decode_debug_.txt");


        for(unsigned int i =0; i<ZigZag_Y_size; i++){
            int val = read_variable_bits();
            ZigZag_Y.push_back(val);
            int rl = read_variable_bits();

            for(int j=0; j<rl; j++){
                ZigZag_Y.push_back(val);
                i++;
            }
            
        }



        for(unsigned int i =0; i<ZigZag_Cb_size; i++){
            int val = read_variable_bits();
            ZigZag_Cb.push_back(val);

            int rl = read_variable_bits();

            for(int j=0; j<rl; j++){
                ZigZag_Cb.push_back(val);
                i++;
            }

        }
        for(unsigned int i =0; i<ZigZag_Cr_size; i++){
            int val = read_variable_bits();
            ZigZag_Cr.push_back(val);

            int rl = read_variable_bits();

            for(int j=0; j<rl; j++){
                ZigZag_Cr.push_back(val);
                i++;
            }

        }
    
        Y = Reverse_ZigZagOrder(ZigZag_Y,height,width);
        Cb_scaled = Reverse_ZigZagOrder(ZigZag_Cb,(height+1)/2,(width+1)/2);
        Cr_scaled = Reverse_ZigZagOrder(ZigZag_Cr,(height+1)/2,(width+1)/2);


        Y = reverse_DCT(Y,height,width,0,q_factor);
        Cb_scaled = reverse_DCT(Cb_scaled,(height+1)/2,(width+1)/2,1,q_factor);
        Cr_scaled = reverse_DCT(Cr_scaled,(height+1)/2,(width+1)/2,1,q_factor);

        auto OG_Y = Y;
        auto OG_Cb = Cb_scaled;
        auto OG_Cr = Cr_scaled;


        if(frame_type){//if p frame

            Y = plus(Y,last_Y);
            Cr_scaled = plus(Cr_scaled,last_Cr);
            Cb_scaled = plus(Cb_scaled,last_Cb);



            for(auto v: motion_vectors){
                int y = v.at(0);
                int x = v.at(1);
                int v_y = v.at(2);
                int v_x = v.at(3);
                
                for(unsigned int k = 0; k<mblock_size; k++){
                    for(unsigned int l = 0; l < mblock_size; l++){
                        Y.at(y+k).at(x+l) = OG_Y.at(y+k).at(x+l) + last_Y.at(y+v_y+k).at(x+v_x+l);
                        Cb_scaled.at((y+k)/2).at((x+l)/2) = OG_Cb.at((y+k)/2).at((x+l)/2) + last_Cb.at(((y+v_y+k)/2)).at(((x+v_x+l)/2));
                        Cr_scaled.at((y+k)/2).at((x+l)/2) = OG_Cr.at((y+k)/2).at((x+l)/2) + last_Cr.at(((y+v_y+k)/2)).at(((x+v_x+l)/2));
                    }
                }
            }
        }

        for (u32 y = 0; y < height; y++)
            for (u32 x = 0; x < width; x++)
                frame.Y(x,y) = clamp(Y.at(y).at(x));
        for (u32 y = 0; y < height/2; y++)
            for (u32 x = 0; x < width/2; x++)
                frame.Cb(x,y) = clamp(Cb_scaled.at(y).at(x));
        for (u32 y = 0; y < height/2; y++)
            for (u32 x = 0; x < width/2; x++)
                frame.Cr(x,y) = clamp(Cr_scaled.at(y).at(x));
        writer.write_frame();

        last_Y = Y;
        last_Cb = Cb_scaled;
        last_Cr = Cr_scaled;
    }


    return 0;
}