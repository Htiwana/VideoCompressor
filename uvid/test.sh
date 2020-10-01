make
time ffmpeg -i ../test_files/flower_cif.y4m -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./uvid_compress 352 288  vlow > flower_compressed.uvi
rm flower_decompressed.y4m
time ./uvid_decompress < flower_compressed.uvi | ffmpeg -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 352x288 -i - -f yuv4mpegpipe flower_decompressed.y4m

time ffmpeg -i ../test_files/harbour_cif.y4m -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./uvid_compress 352 288  low > harbour_compressed.uvi
rm harbour_decompressed.y4m
time ./uvid_decompress < harbour_compressed.uvi | ffmpeg -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 352x288 -i - -f yuv4mpegpipe harbour_decompressed.y4m

time ffmpeg -i ../test_files/news_cif.y4m -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./uvid_compress 352 288  low > news_compressed.uvi
rm news_decompressed.y4m
time ./uvid_decompress < news_compressed.uvi | ffmpeg -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 352x288 -i - -f yuv4mpegpipe news_decompressed.y4m

time ffmpeg -i ../test_files/bunny10.y4m -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./uvid_compress 640 360  low > bunny_compressed.uvi
rm bunny_decompressed.y4m
time ./uvid_decompress < bunny_compressed.uvi | ffmpeg -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 640x360 -i - -f yuv4mpegpipe bunny_decompressed.y4m

time ffmpeg -i ../test_files/elephants10.y4m -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./uvid_compress 640 360  low > elephant_compressed.uvi
rm elephant_decompressed.y4m
time ./uvid_decompress < elephant_compressed.uvi | ffmpeg -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 640x360 -i - -f yuv4mpegpipe elephant_decompressed.y4m

time ffmpeg -i ../test_files/sintel10.y4m -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./uvid_compress 854 480  low > sintel_compressed.uvi
rm sintel_decompressed.y4m
time ./uvid_decompress < sintel_compressed.uvi | ffmpeg -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 854x480 -i - -f yuv4mpegpipe sintel_decompressed.y4m