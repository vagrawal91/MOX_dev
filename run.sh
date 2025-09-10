#matlab -nodesktop -nodisplay -r "cd folder1/; run('mycode.m'); quit"  < /dev/null  > output.txt
#nohup matlab -nodesktop -nodisplay < batch_script_pametricStudy.m > my.log 2>&1 &
nohup /Applications/MATLAB_R2021b.app/bin/matlab -nodesktop -nodisplay < batch_script_pametricStudy.m > my.log 2>&1 &
