Assignment 5 Solution Description

Himmat Singh Tiwana
V00861743


Bit stream format: height u32| width u32| quality 3 bits| frames|

each frame is : number of motion vectors | motion vectors | size of Y channel ( in zig zag flattened form) | Cb size | Cr size | Y channel values | Cb values | Cr values|

everything in frame is stored in variable bits including runs, sizes and actual values. 



	BASIC REQUIREMENTS
	

	I-frames:

my I-frame implementation is mostly a copy of my assignment 4 image compression: 

-An 8x8 block wise DCT for each channel (Y,Cb,Cr)
-Followed by flattening to a ZigZagOrder

- (NEW) after zig zag i have now added rle instead of using delta compression as it seemed to achieve better compression on average ( write_rle() line 231)

- Both the values and run lengths are both written in variable bits and every value is followed by a run length
 ( the stream alternates value, run length, value, run lendth )
 This seemed wasteful at first but because of the ZigZadOrdering and the runs being stored in variable bits small runs are few and stored in few bits.
 
 
	P-frames:
 
I use a frame counter and a simple modulo to switch between I frames and P frames ( line 309). Every I frame is followed by 9 p frames. 

- P blocks or deltas are calculated by subtracting (element wise) the last frame for the current frame
  My P frames effictively have all blocks as delta blocks and no Intra coded blocks. 
  
- the last frame values for each channel are encoded and decoded to have a copy that resembles what the decompressor will receive
	this is done in lines 389 to 397. ( my dct and reverse dct functions include quantization and use q_factor which is determined by the quality setting).
	

	Motion Vectors:

In my pipeline I compute motion vectors after first computing simple deltas for all values in the frame. I am not sure why but this seemed to improve compression performance exponetially. 
For the motion vector compensation I maintain orginal(pre-p-frame-delta-computation) copies of each channel of the current frame. ( OG_Y, OG_Cb, OG_Cr).

motion vectors are computed in find_motion_vectors() line 173. 

I use the fixed neighbour hood search approach. for each 16x16(in Y channel) macroblock I use a fixed search radius of 4. 
If the cumulative squared differences between the pixel values for each channel are lower* for a macroblock in the previous frame only then a motion vector is stored.

*lower than the differences between the current block and the block in the same poisition in the previous frame. ( line 182).

I maintain a vector of motion vectors. Each motion vector is 4 values: the frist two are indices for the block in current frame and the following two are the y axis and x axis deltas to the matching block in the previous frame.(line 197)

The deltas for the blocks with non-zero motion vectors are re-calculated  in motion_vector_deltas() line 211.

	Compression Ratios:
	
I meet and exceed(more in advacned section) the 12x compression ratio requirement on the low quality setting for multiple ( all the ones i trimmed and tried) 10 second clippings from the animated films:

 Elephants Dream (https://media.xiph.org/video/derf/y4m/elephants_dream_360p24.y4m.xz) and 
 Big Buck Bunny https://media.xiph.org/video/derf/y4m/big_buck_bunny_360p24.y4m.xz 
 and the trailer for Sintel https://media.xiph.org/video/derf/y4m/sintel_trailer_2k_480p24.y4m
 as well as harbour_cif.y4m https://media.xiph.org/video/derf/y4m/harbour_cif.y4m
 
  
Since the animated films are too long I used the following commands to get 10 second clippings from different parts of the films:
ffmpeg -ss <30/60/120 etc> -i input.y4m -t 10 output.y4m

	ADVANCED REQUIREMENTS
	
Real Time Decompression: for all the 10 second clips ( harbour, flower, news) as well as 10 secon clippings of the longer 640x360 clips from Elephants Dream and Big Buck Bunny i achieve real time decompression:
all of the clips decode in under 10 seconds on the test server. I used 10 seconds clips of the 640x360 animated clips only to save compression times and i believe the decompression would be real time even for longer snippets of the 640x360 animated clips. 


Significantly Higher Compression ratios: for the animated clips I consistently get 16x compression on average and on some of the 10 second clippings of the Sintel trailer I observe a compression ratio of 36x.

Note: flower and harbour were not quite meeting the 12x compression requirement with my usual q_factor value for low but I was seeing great results of everything else.
So I have added a vlow setting with a higher factor to meet the 12x compression on just those two files with still somewhat resonable quality.


 
 




