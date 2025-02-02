
#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include <tuple>
#include "output_stream.hpp"
#include "yuv_stream.hpp"
#include "uv_common.hpp"
#include <cmath>


OutputBitStream output_stream {std::cout};
const int mblock_size = 16;
const int search_radius =4;

std::vector<std::vector<double>> last_Y;
std::vector<std::vector<double>> last_Cb;
std::vector<std::vector<double>> last_Cr;

std::vector<std::vector<double>> OG_Y;
std::vector<std::vector<double>> OG_Cb;
std::vector<std::vector<double>> OG_Cr;



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

    return (bits*2)+2;
}


//not acutlaly mean just squared diffrences between two macroblocks
double mean_square_diff(int y,int x,int p_y, int p_x){
    double msd=0;

    for(unsigned int k = 0; k<mblock_size; k++){
        for(unsigned int l = 0; l < mblock_size; l++){

            int Y_diff = OG_Y.at(y+k).at(x+l) - last_Y.at(p_y+k).at(p_x+l);
            int Cb_diff = OG_Cb.at((y+k)/2).at((x+l)/2) - last_Cb.at(((p_y+k)/2)).at(((p_x+l)/2));
            int Cr_diff = OG_Cr.at((y+k)/2).at((x+l)/2) - last_Cr.at(((p_y+k)/2)).at(((p_x+l)/2));
            msd+= pow(Y_diff,2);
            msd+= pow(Cb_diff,2);
            msd+= pow(Cr_diff,2);
        }
    }

    return msd;

}

int index_clamp(int index, int dim_max){

    if(index >= (dim_max-mblock_size))
        return dim_max-mblock_size-1;
    else if(index <0 )
        return 0;
    else
        return index;
}

std::vector<std::vector<int>> find_motion_vectors(u32 height, u32 width){

    std::vector<std::vector<int>> motion_vectors;

        //loop over entire frame macro-block by macro-block
        for(int y = 0; y < height-mblock_size; y+=mblock_size){
            for (int x = 0; x < width-mblock_size; x+=mblock_size){

                
                double lowest_sq_diff = mean_square_diff(y,x,y,x);
                std::vector<int> motion_v = {0,0,0,0};

                //for each macroblock check all neighbouring blocks in fixed radius and store motion vec if better than 0
                for(int i = -(mblock_size*search_radius); i<=mblock_size*search_radius; i+=mblock_size){
                    for(int j = -mblock_size*search_radius; j<=mblock_size*search_radius; j+=mblock_size){


                        int v_y = index_clamp(y+i,height);
                        int v_x = index_clamp(x+j,width);
                        
                        double sq_diff = mean_square_diff(y,x,v_y,v_x);

                        if(sq_diff < lowest_sq_diff){
                            lowest_sq_diff = sq_diff;
                            motion_v = { (int)y, (int)x, v_y-y,v_x-x};  
                        }
                    }
                }
                //if non zero motion vector
                if(motion_v.at(2)!=0)
                    motion_vectors.push_back(motion_v);
            }
        }
    
    return motion_vectors;
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> 
motion_vector_deltas(std::vector<std::vector<int>> motion_vectors, std::vector<std::vector<double>> Y,std::vector<std::vector<double>> Cb,std::vector<std::vector<double>> Cr){

    for(auto v: motion_vectors){
        int y = v.at(0);
        int x = v.at(1);
        int v_y = v.at(2);
        int v_x = v.at(3);

        for(unsigned int k = 0; k<mblock_size; k++)
            for(unsigned int l = 0; l < mblock_size; l++){
                Y.at(y+k).at(x+l) = OG_Y.at(y+k).at(x+l) - last_Y.at(y+v_y+k).at(x+v_x+l);
                Cb.at((y+k)/2).at((x+l)/2) = OG_Cb.at((y+k)/2).at((x+l)/2) - last_Cb.at(((y+v_y+k)/2)).at(((x+v_x+l)/2));
                Cr.at((y+k)/2).at((x+l)/2) = OG_Cr.at((y+k)/2).at((x+l)/2) - last_Cr.at(((y+v_y+k)/2)).at(((x+v_x+l)/2));
            }
                
    }

    return {Y,Cb,Cr};
}

int write_rle(std::vector<double> vec){

    int rl =0;
    write_variable_bits(vec.at(0));
    for(unsigned int i =1; i<vec.size(); i++){


        if(vec.at(i) == vec.at(i-1)){
            rl++;
        }else{
            write_variable_bits(rl);
            write_variable_bits(vec.at(i));
            rl=0;
        }
        
    }
    write_variable_bits(rl);

    return 0;
}



int main(int argc, char** argv){

    std::ofstream debugfile;
    debugfile.open ("DEBUG.txt");

    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <width> <height> <vlow/low/medium/high>" << std::endl;
        return 1;
    }

    u32 width = std::stoi(argv[1]);
    u32 height = std::stoi(argv[2]);
    std::string quality{argv[3]};

    YUVStreamReader reader {std::cin, width, height};


    output_stream.push_u32(height);
    output_stream.push_u32(width);


    unsigned int quantization_factor = 0;

    if(quality=="vlow"){
        quantization_factor = 4;
    }else if(quality =="low"){
        quantization_factor =1;
    }else if(quality == "medium"){
        quantization_factor =2;
    }else{
        quantization_factor =3;
    }
    
    output_stream.push_bits(quantization_factor,3);
    
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
        case 4:
            q_factor = 6;
            break;
    }

    unsigned int frame_type =0;
    unsigned int frame_num =0;

    last_Y = create_2d_vector<double>(height,width);
    last_Cb = create_2d_vector<double>(height/2,width/2);
    last_Cr = create_2d_vector<double>(height/2,width/2);

    while (reader.read_next_frame()){

        frame_num++;
        frame_type = (frame_num%10==0)?0:1;// 0 = I frame, 1 = P frame

        output_stream.push_byte(1); //Use a one byte flag to indicate whether there is a frame here
        
        YUVFrame420& frame = reader.frame();
        output_stream.push_bits(frame_type,2);
        
        auto Y = create_2d_vector<double>(height,width);
        auto Cb = create_2d_vector<double>(height/2,width/2);
        auto Cr = create_2d_vector<double>(height/2,width/2);

        //Extract the planes into different arrays 
        OG_Y = create_2d_vector<double>(height,width);
        OG_Cb = create_2d_vector<double>(height/2,width/2);
        OG_Cr = create_2d_vector<double>(height/2,width/2);
        for(unsigned int y = 0; y < height; y++)
            for (unsigned int x = 0; x < width; x++)
                OG_Y.at(y).at(x) = frame.Y(x,y);
        for(unsigned int y = 0; y < height/2; y++)
            for (unsigned int x = 0; x < width/2; x++)
                OG_Cb.at(y).at(x) = frame.Cb(x,y);
        for(unsigned int y = 0; y < height/2; y++)
            for (unsigned int x = 0; x < width/2; x++)
                OG_Cr.at(y).at(x) = frame.Cr(x,y);


        

        auto motion_vectors = find_motion_vectors(height,width);


        if(frame_type){//if p frame

            //lines 361 to 363 seem to fix accumulating artifacts that look like delta value deterioration ....seems excessive to encode/decode current frame too but it works
            OG_Y = reverse_DCT(DCT(OG_Y,height,width,0,q_factor),height,width,0,q_factor);
            OG_Cb = reverse_DCT(DCT(OG_Cb,height/2,width/2,0,q_factor),height/2,width/2,0,q_factor);
            OG_Cr = reverse_DCT(DCT(OG_Cr,height/2,width/2,0,q_factor),height/2,width/2,0,q_factor);  

            //compute deltas for all values in frame ( same effect as doing all blocks, block by block)
            Y = minus(OG_Y,last_Y);
            Cb = minus(OG_Cb,last_Cb);
            Cr = minus(OG_Cr,last_Cr);

            //compute macro-block wise motion vectors 
            auto motion_vectors = find_motion_vectors(height,width);

            //recompute deltas for macro-blocks with better than 0 motion vectors
            std::tie(Y,Cb,Cr) = motion_vector_deltas(motion_vectors,Y,Cb,Cr);

            //write motion vectors
            write_variable_bits(motion_vectors.size());
            for(unsigned int i =0; i<motion_vectors.size(); i++)
                for(unsigned int j =0; j<4; j++)
                    write_variable_bits(motion_vectors.at(i).at(j));

        }else{
            Y = OG_Y;
            Cb = OG_Cb;
            Cr = OG_Cr;
        }

     
        //DCT for each plane
        auto DCT_Y = DCT(Y,height,width,0, q_factor);     
        auto DCT_Cb = DCT(Cb,(height+1)/2,(width+1)/2,1,q_factor);
        auto DCT_Cr = DCT(Cr,(height+1)/2,(width+1)/2,1,q_factor);

        //flatten each channel in zig zag order
        auto ZigZag_Y = ZigZagOrder(DCT_Y,height,width);        
        auto ZigZag_Cb = ZigZagOrder(DCT_Cb,(height+1)/2,(width+1)/2);
        auto ZigZag_Cr = ZigZagOrder(DCT_Cr,(height+1)/2,(width+1)/2);

        write_variable_bits(ZigZag_Y.size());
        write_variable_bits(ZigZag_Cb.size());
        write_variable_bits(ZigZag_Cr.size());

        //Write ZigZag values in variable bit format with rle
        write_rle(ZigZag_Y);
        write_rle(ZigZag_Cb);
        write_rle(ZigZag_Cr);

        //compute last frame vlaues that decompressor will obtain by doing and undoing DCT on orginal values
        OG_Y = DCT(OG_Y,height,width,0,q_factor);
        last_Y = reverse_DCT(OG_Y,height,width,0,q_factor);

        OG_Cb = DCT(OG_Cb,(height+1)/2,(width+1)/2,1,q_factor);
        last_Cb = reverse_DCT(OG_Cb,(height+1)/2,(width+1)/2,1,q_factor);

        OG_Cr = DCT(OG_Cr,(height+1)/2,(width+1)/2,1,q_factor);
        last_Cr = reverse_DCT(OG_Cr,(height+1)/2,(width+1)/2,1,q_factor);
    }

    output_stream.push_byte(0); //Flag to indicate end of data
    output_stream.flush_to_byte();
    debugfile.close();

    return 0;
}