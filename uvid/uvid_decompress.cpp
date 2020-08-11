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
#include "uvg_common.hpp"


InputBitStream input_stream {std::cin};

std::vector<std::vector<double>> DCT_matrix  = {{ 0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536}, 
                                                { 0.4904,  0.4157,  0.2778,  0.0975, -0.0975, -0.2778, -0.4157, -0.4904}, 
                                                { 0.4619,  0.1913, -0.1913, -0.4619, -0.4619, -0.1913,  0.1913,  0.4619}, 
                                                { 0.4157, -0.0975, -0.4904, -0.2778,  0.2778,  0.4904,  0.0975, -0.4157},
                                                { 0.3536, -0.3536, -0.3536,  0.3536,  0.3536, -0.3536, -0.3536,  0.3536},
                                                { 0.2778, -0.4904,  0.0975,  0.4157, -0.4157, -0.0975,  0.4904, -0.2778},
                                                { 0.1913, -0.4619,  0.4619, -0.1913, -0.1913,  0.4619, -0.4619,  0.1913},
                                                { 0.0975, -0.2778,  0.4157, -0.4904,  0.4904, -0.4157,  0.2778, -0.0975}};

std::vector<std::vector<double>> Y_quant  = {{ 16, 11, 10, 16, 24, 40, 51, 61}, 
                            { 12, 12, 14, 19, 26, 58, 60, 55}, 
                            { 14, 13, 16, 24, 40, 57, 69, 56}, 
                            { 14, 17, 22, 29, 51, 87, 80, 62},
                            { 18, 22, 37, 56, 68, 109, 103, 77},
                            { 24, 35, 55, 64, 81, 104, 113, 92},
                            { 49, 64, 78, 87, 103, 121, 120, 101},
                            { 72, 92, 95, 98, 112, 100, 103, 99}};

std::vector<std::vector<double>> C_quant  ={{17, 18, 24, 47, 99, 99, 99, 99},
                        {18, 21, 26, 66, 99, 99, 99, 99},
                        {24, 26, 56, 99, 99, 99, 99, 99},
                        {47, 66, 99, 99, 99, 99, 99, 99},
                        {99, 99, 99, 99, 99, 99, 99, 99},
                        {99, 99, 99, 99, 99, 99, 99, 99},
                        {99, 99, 99, 99, 99, 99, 99, 99},
                        {99, 99, 99, 99, 99, 99, 99, 99}};

//matrix multiplication from https://www.programiz.com/cpp-programming/examples/matrix-multiplication
std::vector<std::vector<double>> multiply( std::vector<std::vector<double>> a, std::vector<std::vector<double>> b){

    auto mult = create_2d_vector<double>(8,8);

    for(int i = 0; i < 8; ++i)
        for(int j = 0; j < 8; ++j)
            for(int k = 0; k < 8; ++k)
            {
                mult[i][j] += a[i][k] * b[k][j];
            }

    return mult;
}

std::vector<std::vector<double>> transpose( std::vector<std::vector<double>> a){

    auto transpose = create_2d_vector<double>(8,8);

    for(int i = 0; i < 8; ++i)
        for(int j = 0; j < 8; ++j) {
            transpose[j][i] = a[i][j];
        }

    return transpose;
}

std::vector<std::vector<double>> reverse_DCT(std::vector<std::vector<double>> plane, unsigned int height, unsigned int width, int channel, double q_factor){
    std::vector<std::vector<double>> quantizer =  ((channel==0)?Y_quant:C_quant);

    auto block = create_2d_vector<double>(8,8);
    std::vector<std::vector<double>> transformation = create_2d_vector<double>(height,width);

    //read as 8x8 blocks
    double val = 0;
    for(unsigned int y = 0; y < height; y+=8){
        for (unsigned int x = 0; x < width; x+=8){
            for(unsigned int i = 0; i < 8; i++){
                for(unsigned int j = 0; j < 8; j++){ 
                    if(((y+i)>=height || (x+j) >=width))//only updates within image and previous pixel value gets reused outside. so takes care of padding ? 
                        break;
                    val = plane.at(y+i).at(x+j);   
                    block.at(i).at(j) = val*(quantizer[i][j]*q_factor);
                }
            }
            //DCT on block
            auto DCT = multiply(transpose(DCT_matrix),block);
            DCT = multiply(DCT,DCT_matrix);
            for(unsigned int i = 0; i < 8; i++){
                for(unsigned int j = 0; j < 8; j++){
                    if((y+i)>=height || (x+j) >= width)
                        break; 
                    transformation.at(y+i).at(x+j) = std::round(DCT.at(i).at(j));
                }
            }
        }
    }

    return transformation;
}


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

unsigned char clamp (double val){

    if(val > 255){
        return 255;
    }else if(val < 0){
        return 0;
    }else{
        return (unsigned char)(unsigned int)val;
    }
}

std::vector<std::vector<double>> plus(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B){
    std::vector<std::vector<double>> result = create_2d_vector<double>(A.size(),A.at(0).size());

    for(unsigned int i =0; i<A.size(); i++)
        for(unsigned int j =0; j<A.at(0).size(); j++)
            result.at(i).at(j) = A.at(i).at(j) + B.at(i).at(j);

    return result;
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

    //Note: This program must not take any command line arguments. (Anything
    //      it needs to know about the data must be encoded into the bitstream)
    
    

    u32 height {input_stream.read_u32()};
    u32 width {input_stream.read_u32()};

    YUVStreamWriter writer {std::cout, width, height};

    auto last_Y = create_2d_vector<double>(height,width);
    auto last_Cb = create_2d_vector<double>(height/2,width/2);
    auto last_Cr = create_2d_vector<double>(height/2,width/2);
    
    unsigned int frame_num =0;

    while (input_stream.read_byte()){
        YUVFrame420& frame = writer.frame();
        frame_num++;

        unsigned int frame_type = input_stream.read_bits(2);

        
        //motion vectors
        unsigned int motion_vector_size = read_variable_bits();
        std::vector<std::vector<int>> motion_vectors;
        std::vector<int> motion_v;

        for(int i =0; i<motion_vector_size; i++){
            for(int j =0; j<4; j++){
                motion_v.push_back(read_variable_bits());
            }
            motion_vectors.push_back(motion_v);
            motion_v.clear();
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

            //rle
            // if(val==0){
            //     int rl = read_variable_bits();
            //     if(frame_num==2)
            //         debugfile << "run: " << rl << std::endl;
            //     for(int j =0; j<rl; j++){
            //         ZigZag_Y.push_back(val);
            //         i++;
            //     }
                
            // }

            
        }

        // for(unsigned int i =1; i<ZigZag_Y.size(); i++){
        //     if(i%64>6)
        //         ZigZag_Y.at(i)+=ZigZag_Y.at(i-1);
        // }

        for(unsigned int i =0; i<ZigZag_Cb_size; i++){
            int val = read_variable_bits();
            ZigZag_Cb.push_back(val);
        }
        for(unsigned int i =0; i<ZigZag_Cr_size; i++){
            int val = read_variable_bits();
            ZigZag_Cr.push_back(val);
        }
    
        Y = Reverse_ZigZagOrder(ZigZag_Y,height,width);
        Cb_scaled = Reverse_ZigZagOrder(ZigZag_Cb,(height+1)/2,(width+1)/2);
        Cr_scaled = Reverse_ZigZagOrder(ZigZag_Cr,(height+1)/2,(width+1)/2);


        Y = reverse_DCT(Y,height,width,0,1);
        Cb_scaled = reverse_DCT(Cb_scaled,(height+1)/2,(width+1)/2,1,1);
        Cr_scaled = reverse_DCT(Cr_scaled,(height+1)/2,(width+1)/2,1,1);

        auto OG_Y = Y;
        if(frame_type){//if p frame
            Y = plus(Y,last_Y);
            Cr_scaled = plus(Cr_scaled,last_Cr);
            Cb_scaled = plus(Cb_scaled,last_Cb);
        }

        for(auto v: motion_vectors){
            int y = v.at(0);
            int x = v.at(1);
            int v_y = v.at(2);
            int v_x = v.at(3);
            
            for(unsigned int k = 0; k<8; k++){
                for(unsigned int l = 0; l < 8; l++){
                    Y.at(y+k).at(x+l) = OG_Y.at(y+k).at(x+l)+last_Y.at(y+v_y+k).at(x+v_x+l);
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