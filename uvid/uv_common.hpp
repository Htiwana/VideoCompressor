/* uvg_common.hpp

   B. Bird - 07/02/2020

   modified by Himmat Singh Tiwana 08/12/2020


*/

#ifndef UV_COMMON_HPP
#define UV_COMMON_HPP

#include <vector>
#include <string>


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


//Convenience function to wrap around the nasty notation for 2d vectors
template<typename T>
std::vector<std::vector<T> > create_2d_vector(unsigned int outer, unsigned int inner){
    std::vector<std::vector<T> > V {outer, std::vector<T>(inner,T() )};
    return V;
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
                    transformation.at(y+i).at(x+j) = DCT.at(i).at(j);
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

unsigned char clamp (double val){

    if(val > 255){
        return 255;
    }else if(val < 0){
        return 0;
    }else{
        return (unsigned char)(unsigned int)val;
    }
}



//The floating point calculations we use while converting between 
//RGB and YCbCr can occasionally yield values slightly out of range
//for an unsigned char (e.g. -1 or 255.9).
//Furthermore, we want to ensure that any conversion uses rounding
//and not truncation (to improve accuracy).
inline unsigned char round_and_clamp_to_char(double v){
    //Round to int 
    int i = (int)(v+0.5);
    //Clamp to the range [0,255]
    if (i < 0)
        return 0;
    else if (i > 255)
        return 255;
    return i;
}

/* The exact RGB <-> YCbCr conversion formula used here is the "JPEG style"
   conversion (there is some debate over the best conversion formula) */
struct PixelYCbCr;
struct PixelRGB{
    unsigned char r, g, b;
    PixelYCbCr to_ycbcr(); //Implementation is below (since the PixelYCbCr type has to exist before we can fully define this function)
};

struct PixelYCbCr{
    unsigned char Y, Cb, Cr;
    inline PixelRGB to_rgb(){
        return {
            round_and_clamp_to_char(Y + 1.402*(Cr-128.0)),
            round_and_clamp_to_char(Y-0.344136*(Cb-128.0)-0.714136*(Cr-128.0)),
            round_and_clamp_to_char(Y+1.772*(Cb-128.0))
        };
    }
};


inline PixelYCbCr PixelRGB::to_ycbcr(){
    return {
        round_and_clamp_to_char(0.299*r + 0.587*g + 0.114*b),
        round_and_clamp_to_char(128 + -0.168736*r + -0.331264*g + 0.5*b),
        round_and_clamp_to_char(128 + 0.5*r + -0.418688*g + -0.081312*b)
    };
}


#endif