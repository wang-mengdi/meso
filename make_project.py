import os
import sys
import platform

if __name__=='__main__':
    assert (len(sys.argv)==2), 'Error: must pass 1 parameter(s). Example: python make_project.py _tests_reservoir'
    proj_name=sys.argv[1]
    bin_dir=os.path.join('..','..','bin',proj_name)
    build_dir=os.path.join('..','..','build',proj_name)
    lua_file=os.path.join('proj',proj_name,'xmake.lua')
    proj_dir=os.path.join('proj',proj_name)
    clean_cmd="xmake c -P {} -a".format(proj_dir)
    config_cmd="xmake f -o {} -P {}".format(bin_dir, proj_dir)
    project_cmd="xmake project -k vsxmake -v -a \"x64\" -P {} {}".format(proj_dir, build_dir)
    print(clean_cmd)
    os.system(clean_cmd)
    print(config_cmd)
    os.system(config_cmd)
    print(project_cmd)
    os.system(project_cmd)