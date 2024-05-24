import os
import sys
import platform

if __name__=='__main__':
    proj_name="openvdb-test"
    if len(sys.argv)==2:
        proj_name=sys.argv[1]
    #bin_dir=os.path.join('..','..','bin',proj_name)
    bin_dir=os.path.join('bin',proj_name)
    lua_file=os.path.join('proj',proj_name,'xmake.lua')
    proj_dir=os.path.join('proj',proj_name)
    clean_cmd="xmake c -P {} -a".format(proj_dir)
    config_cmd="xmake f -o {} -P {} -y".format(bin_dir, proj_dir)
    build_cmd = "xmake build -v -P {} {}".format(proj_dir, proj_name)
    print(clean_cmd)
    os.system(clean_cmd)
    print(config_cmd)
    os.system(config_cmd)
    print(build_cmd)
    os.system(build_cmd)