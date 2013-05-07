Install boost v1.53:
On MAX OS:
install_name_tool -id /Users/pengh/Documents/boost/lib/libboost_system.dylib libboost_system.dylib
install_name_tool -change libboost_system.dylib @executable_path/../Frameworks/libboost_system.dylib libboost_filesystem.dylib
For old versions:
remove -lboost_system

Install samtools v0.1.19:
