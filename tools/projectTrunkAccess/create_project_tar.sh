#!/bin/bash
#Creates a compressed project tarball containing a version of FLASH trunk
#that can be distributed.  All directories that appear in trunk_additions
#file are removed unless they appear in trunk_access file.  The directories
#listed in trunk_access file are specific to a particular project, e.g.
#Northwestern folks need access to the new I/O code in FLASH trunk.

#Parameters.
#-----------------------------------------------------------------------
flash_trunk="flash.uchicago.edu/home/svn/repos/FLASH3/trunk"
revision=13506
user="cdaley"
check_out_path="/home/cdaley/flash"

check_out_name="FLASH3_Mats_r${revision}"
#check_out_name="FLASH3_Northwestern_r${revision}"

trunk_additions="trunk_additions.txt"

trunk_access="trunk_access_mats.txt"
#trunk_access="trunk_access_northwestern.txt"


#Do not edit below this line.
#-----------------------------------------------------------------------
start_dir=$(pwd)
work_dir="${check_out_path}/${check_out_name}"
source_tree_info="source_tree_info.txt"
source_tree_deletions="source_tree_deletions.txt"
tmp_trunk_access="sorted_${trunk_access}"

echo "Checking out r${revision} of ${flash_trunk} as ${work_dir}"
svn co -r "${revision}" svn+ssh://"${user}@${flash_trunk}" "${work_dir}"
if [ $? -ne 0 ]; then
    exit 1
fi

cd "${work_dir}"
#Use $(pwd) in echo statement just in case we are not in work_dir.
echo -e "Entered directory $(pwd)/\nStart deleting from this directory? [y|n]"
read user_continue

if [ "${user_continue}" == "y" ]; then

    echo -e "Creating pruned FLASH source tree" 2>&1 | tee "${source_tree_info}"

    #Construct a list of directories that should be deleted from FLASH trunk.
    #Use a temporary file to ensure data is sorted for "comm" command.
    cat "${start_dir}/${trunk_access}" | sort > "${tmp_trunk_access}"
    comm -13 "${tmp_trunk_access}" "${start_dir}/${trunk_additions}" \
	> "${source_tree_deletions}"
    rm -f "${tmp_trunk_access}"

    #Delete docs to be on the safe side.
    echo "docs" >> "${source_tree_deletions}"


    #Delete the list of directories from the checked out copy of FLASH trunk.
    while read line; do
	echo "Removing ${line}" 2>&1 | tee -a "${source_tree_info}"
	rm -rf "${line}"
    done <"${source_tree_deletions}"
    rm -f "${source_tree_deletions}"


    #Svn information -- Use '!' to omit information about removed code.
    echo -e "\nExtracting various svn information before .svn removal" \
	2>&1 | tee -a "${source_tree_info}"
    svn info >> "${source_tree_info}"
    svn st -v | grep -Fv ! >> "${source_tree_info}"

    #The -depth argument in this long command ensures parent directories are
    #not deleted before child directories.  For some reason it is only needed
    #in bash scripts and not on the command line?
    find . -depth -name '.svn' -type d -exec rm -rf {} \;


    echo -e "\nCreating compressed tar archive named ${check_out_name}.tar.gz"
    cd ..
    tar -cf "${check_out_name}".tar "${check_out_name}/"
    gzip -9 "${check_out_name}".tar
    rm -rf "${check_out_name}/"
fi

echo -e "\nReturning from $(pwd) to ${start_dir}"
cd "${start_dir}"
echo -e "\nEnd of script" 
