#!/bin/bash

while getopts ":r:i:o:I:h:n:" opt; do
  case $opt in
    
    h) 
    echo "       Help"
    echo ""
    echo "       Script runs docker base container for NGS pipelines."
    echo ""
    echo "      ************ OPTIONS ****************"
    echo "      -r: Reference dir "
    echo "      -i: Input docker mount point folder"
    echo "      -o: Output docker mount point folder"
    echo "      -I: Base container image"
    echo "      -n: Container name (username)*"
    echo "       *make sure that container name is not being used,"
    echo "        try docker ps -a |grep [container name]"
    echo "      *************************************" 
	echo ""
    echo "       COMMAND SYNTAX:"
    echo "       (all options MUST be specified)"
    echo ""
    echo "       bash run_ngs_container.sh -r [reference dir] -i [data input dir] -o [data output dir] -I [base container image]"
    echo""
    echo "       End of help."
	echo ""
    exit
      ;;
    r) reference_folder="$OPTARG"	
    echo "Refernce folder: $OPTARG" >&2
      ;;    
    i) input_folder="$OPTARG"
    echo "Input docker mount point folder: $OPTARG" >&2
      ;;
    o) output_folder="$OPTARG"
    echo "Output docker mount point folder: $OPTARG" >&2
      ;;
    I) docker_image="$OPTARG"
    echo "Base container docker image: $OPTARG" >&2
      ;;  
   
    n) container_name="$OPTARG"
    echo "Container name: $OPTARG" >&2
      ;;      
    \?)
      echo "Invalid option: $OPTARG , run -h help " >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument. Run -h help" >&2
      exit 1
      ;;
  esac
done

echo ""
echo "starting docker container..."
docker run --name $container_name -v /var/run/docker.sock:/var/run/docker.sock -v $input_folder/:$input_folder/ -v $output_folder://media/DISK1/exome_aggregation/exome_test/OUTPUT/ -v $reference_folder://media/DISK1/exome_aggregation/GATK_b37 -it $docker_image

