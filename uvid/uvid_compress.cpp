/* uvid_compress.cpp
   CSC 485B/578B/SENG 480B - Data Compression - Summer 2020

   Starter code for Assignment 5

   Reads video data from stdin in uncompresed YCbCr (YUV) format 
   (With 4:2:0 subsampling). To produce this format from 
   arbitrary video data in a popular format, use the ffmpeg
   tool and a command like 

     ffmpeg -i videofile.mp4 -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./this_program <width> height>

   Note that since the width/height of each frame is not encoded into the raw
   video stream, those values must be provided to the program as arguments.

   B. Bird - 07/15/2020
*/

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include <tuple>
#include "output_stream.hpp"
#include "yuv_stream.hpp"
#include "bitmap_image.hpp"
#include "uvg_common.hpp"

OutputBitStream output_stream {std::cout};

//https://www.mathworks.com/help/images/ref/dctmtx.html
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



//my zig zag implementation was influenced by https://www.geeksforgeeks.org/zigzag-or-diagonal-traversal-of-matrix/
//I did not actually understand or copy their implementation exactly but mainly borrowed the idea of doing it in two sets of nested loops ( once going along the fist column and then last row)
//instead of trying to do it in one set of nested loops
std::vector<double> ZigZagBlock(std::vector<std::vector<double>> block){

    std::vector<double> ZigZag;

    int row =0;
    int column =0;
    int inc = 1;

    //zig zag till half ( longest diagonal starting from index 7,0)
    for(unsigned int i =0; i<block.size(); i++){
        row =(i%2==0)?0:i;
        column = (i%2==0)?i:0;
        inc =(i%2==0)?-1:1;

        ZigZag.push_back(block.at(row).at(column));
        for(unsigned int j =0; j<i; j++){
            row+=-inc;
            column+=inc;
            ZigZag.push_back(block.at(row).at(column));
        }
    }

    //zig from longest dianoal to end ( index 7,7)
    for(unsigned int i =1; i<block.size(); i++){
        row =(i%2==0)?7:i;
        column = (i%2==0)?i:7;
        inc =(i%2==0)?1:-1;


        for(unsigned int j =i; j<block.size(); j++){
            ZigZag.push_back(block.at(row).at(column));
            row+=-inc;
            column+=inc;
        }
        
    }

    return ZigZag;


}

std::vector<double> ZigZagOrder(std::vector<std::vector<double>> plane, unsigned int height, unsigned int width){


    auto block = create_2d_vector<double>(8,8);
    std::vector<double> ZigZag;

    //read as 8x8 blocks
    double val = 0;
    for(unsigned int y = 0; y < height; y+=8){
        for (unsigned int x = 0; x < width; x+=8){
            for(unsigned int i = 0; i < 8; i++){
                for(unsigned int j = 0; j < 8; j++){
                    if(!((y+i)>=height || (x+j) >=width))
                        val = plane.at(y+i).at(x+j);
                    block.at(i).at(j) = val;
                }
            }
            //turn each block into zig zag order and append
            auto ZigZag_block = ZigZagBlock(block);
            
            ZigZag.insert(ZigZag.end(),ZigZag_block.begin(),ZigZag_block.end());
        }
    }

    return ZigZag;

}


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


//reads plane in 8x8 blocks,  performs DCT on each blocks and puts quantized values in transformation
std::vector<std::vector<double>> DCT(std::vector<std::vector<double>> plane, unsigned int height, unsigned int width, int channel, double q_factor){
    std::vector<std::vector<double>> quantizer =  ((channel==0)?Y_quant:C_quant);

    auto block = create_2d_vector<double>(8,8);
    std::vector<std::vector<double>> transformation = create_2d_vector<double>(height,width);

    //read as 8x8 blocks
    double val = 0;
    for(unsigned int y = 0; y < height; y+=8){
        for (unsigned int x = 0; x < width; x+=8){
            for(unsigned int i = 0; i < 8; i++){
                for(unsigned int j = 0; j < 8; j++){ 
                    if(!((y+i)>=height || (x+j) >=width))//only updates within image and previous pixel value gets reused outside. so takes care of padding ? 
                        val = plane.at(y+i).at(x+j);   
                    block.at(i).at(j) = val;
                }
            }
            //DCT on block
            auto DCT = multiply(DCT_matrix,block);
            DCT = multiply(DCT,transpose(DCT_matrix));
            for(unsigned int i = 0; i < 8; i++){
                for(unsigned int j = 0; j < 8; j++){
                    if((y+i)>=height || (x+j) >= width)
                        break; 
                    transformation.at(y+i).at(x+j) = std::round(DCT.at(i).at(j)/(quantizer[i][j]*q_factor));
                }
            }
        }
    }

    return transformation;
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

std::vector<std::vector<double>> minus(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B){
    std::vector<std::vector<double>> result = create_2d_vector<double>(A.size(),A.at(0).size());

    assert(A.size()==B.size());
    for(unsigned int i =0; i<A.size(); i++)
        for(unsigned int j =0; j<A.at(0).size(); j++)
            result.at(i).at(j) = A.at(i).at(j) - B.at(i).at(j);

    return result;
}




std::vector<std::vector<double>> plus(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B){
    std::vector<std::vector<double>> result = create_2d_vector<double>(A.size(),A.at(0).size());

    for(unsigned int i =0; i<A.size(); i++)
        for(unsigned int j =0; j<A.at(0).size(); j++)
            result.at(i).at(j) = A.at(i).at(j) + B.at(i).at(j);

    return result;
}


int write_variable_bits(int val){

    unsigned int sign = (val<0)?1:0;

    output_stream.push_bit(sign);
    val = (val<0)?(val*-1):val;

    unsigned int bits = ceil(log2(val));
    bits = (!(val & (val -1)))?bits+1:bits;//if power of two then one more bit than log 2 https://stackoverflow.com/questions/108318/whats-the-simplest-way-to-test-whether-a-number-is-a-power-of-2-in-c
    
    for(unsigned int i = 0; i < bits; i++)
        output_stream.push_bit(1);
    output_stream.push_bit(0);
    output_stream.push_bits((unsigned int)val,bits);

    return bits;
}


int rle(int index, std::vector<double> vec){

    int cur_block = index/64;
    int block_end = (cur_block+1)*64;

    int run_length=0;
    for(int i = index+1; i< block_end; i++){
        if(vec.at(i)==vec.at(index)){
            run_length++;
        }else{
            return run_length;
        }
    }

    return run_length;
}

int main(int argc, char** argv){

    std::ofstream debugfile;
    debugfile.open ("debug_.txt");

    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <width> <height> <low/medium/high>" << std::endl;
        return 1;
    }
    u32 width = std::stoi(argv[1]);
    u32 height = std::stoi(argv[2]);
    std::string quality{argv[3]};

    YUVStreamReader reader {std::cin, width, height};


    output_stream.push_u32(height);
    output_stream.push_u32(width);

    unsigned int frame_type =0;
    unsigned int frame_num =0;

    auto last_Y = create_2d_vector<double>(height,width);
    auto last_Cb = create_2d_vector<double>(height/2,width/2);
    auto last_Cr = create_2d_vector<double>(height/2,width/2);

    while (reader.read_next_frame()){

        frame_num++;
        frame_type = (frame_num%5==0)?0:1;// 0 = I frame, 1 = P frame
        //frame_type =0;

        output_stream.push_byte(1); //Use a one byte flag to indicate whether there is a frame here
        
        YUVFrame420& frame = reader.frame();
        output_stream.push_bits(frame_type,2);

        
        auto Y = create_2d_vector<double>(height,width);
        auto Cb = create_2d_vector<double>(height/2,width/2);
        auto Cr = create_2d_vector<double>(height/2,width/2);

        //Extract the Y plane into its own array 
        auto OG_Y = create_2d_vector<double>(height,width);
        for(unsigned int y = 0; y < height; y++)
            for (unsigned int x = 0; x < width; x++)
                OG_Y.at(y).at(x) = frame.Y(x,y);
        
        //Extract the Cb plane into its own array 
        auto OG_Cb = create_2d_vector<double>(height/2,width/2);
        for(unsigned int y = 0; y < height/2; y++)
            for (unsigned int x = 0; x < width/2; x++)
                OG_Cb.at(y).at(x) = frame.Cb(x,y);

        //Extract the Cr plane into its own array
        auto OG_Cr = create_2d_vector<double>(height/2,width/2);
        for(unsigned int y = 0; y < height/2; y++)
            for (unsigned int x = 0; x < width/2; x++)
                OG_Cr.at(y).at(x) = frame.Cr(x,y);

        
        if(frame_type){//if p frame        
            Y = minus(OG_Y,last_Y);
            Cb = minus(OG_Cb,last_Cb);
            Cr = minus(OG_Cr,last_Cr);
        }else{
            Y = OG_Y;
            Cb = OG_Cb;
            Cr = OG_Cr;
        }
        
        //DCT for each plane
        auto DCT_Y = DCT(Y,height,width,0,1);     
        auto DCT_Cb = DCT(Cb,(height+1)/2,(width+1)/2,1,1);
        auto DCT_Cr = DCT(Cr,(height+1)/2,(width+1)/2,1,1);


        auto ZigZag_Y = ZigZagOrder(DCT_Y,height,width);        
        auto ZigZag_Cb = ZigZagOrder(DCT_Cb,(height+1)/2,(width+1)/2);
        auto ZigZag_Cr = ZigZagOrder(DCT_Cr,(height+1)/2,(width+1)/2);

        write_variable_bits(ZigZag_Y.size());
        write_variable_bits(ZigZag_Cb.size());
        write_variable_bits(ZigZag_Cr.size());


        //compute deltas
        // for(unsigned int i =1; i<ZigZag_Y.size(); i++){
        //     if(i%64>4)
        //         ZigZag_Y.at(i) = ZigZag_Y.at(i)-ZigZag_Y.at(i-1);
        // }

        //Write ZigZag values in variable bit format
        for(unsigned int i =0; i<ZigZag_Y.size(); i++){
            int val = std::round(ZigZag_Y.at(i));
            write_variable_bits(val);

            //rle
            // if(val==0){
            //     int rl = rle(i,ZigZag_Y);
            //     if(frame_num==2)
            //         debugfile2 << "run: " << rl << std::endl;
            //     write_variable_bits(rl);
            //     i+= rl;
            // }
        }

        for(unsigned int i =0; i<ZigZag_Cb.size(); i++)
            write_variable_bits(std::round(ZigZag_Cb.at(i)));
        for(unsigned int i =0; i<ZigZag_Cr.size(); i++)
            write_variable_bits(std::round(ZigZag_Cr.at(i)));

        //compute last frame vlaues that decompressor will see by doing and undoing DCT
        OG_Y = DCT(OG_Y,height,width,0,1);
        last_Y = reverse_DCT(OG_Y,height,width,0,1);

        OG_Cb = DCT(OG_Cb,(height+1)/2,(width+1)/2,1,1);
        last_Cb = reverse_DCT(OG_Cb,(height+1)/2,(width+1)/2,1,1);

        OG_Cr = DCT(OG_Cr,(height+1)/2,(width+1)/2,1,1);
        last_Cr = reverse_DCT(OG_Cr,(height+1)/2,(width+1)/2,1,1);

    }

    output_stream.push_byte(0); //Flag to indicate end of data
    output_stream.flush_to_byte();
    debugfile.close();

    return 0;
}