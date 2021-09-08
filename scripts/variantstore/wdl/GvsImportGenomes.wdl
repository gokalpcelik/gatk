version 1.0

workflow GvsImportGenomes {

  input {
    Array[File] input_vcfs
    Array[String] external_sample_names
    File interval_list
    String output_directory
    File? sample_map
    String project_id
    String dataset_name
    String? pet_schema = "location:INTEGER,sample_id:INTEGER,state:STRING"
    String? vet_schema = "sample_id:INTEGER,location:INTEGER,ref:STRING,alt:STRING,AS_RAW_MQ:STRING,AS_RAW_MQRankSum:STRING,QUALapprox:STRING,AS_QUALapprox:STRING,AS_RAW_ReadPosRankSum:STRING,AS_SB_TABLE:STRING,AS_VarDP:STRING,call_GT:STRING,call_AD:STRING,call_GQ:INTEGER,call_PGT:STRING,call_PID:STRING,call_PL:STRING"
    String? service_account_json_path
    String? drop_state = "SIXTY"
    Boolean? drop_state_includes_greater_than = false

    Int batch_size = 1

    Int? preemptible_tries
    File? gatk_override
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])

  # return an error if the lengths are not equal
  Int input_length = length(input_vcfs)
  if (input_length != length(external_sample_names)) {
    call TerminateWorkflow {
      input:
        message = 'Input array lengths do not match'
    }
  }

  call SetLock {
    input:
      output_directory = output_directory,
      service_account_json_path = service_account_json_path,
      preemptible_tries = preemptible_tries
  }

  call GetSampleIds {
    input:
      external_sample_names = external_sample_names,
      project_id = project_id,
      dataset_name = dataset_name,
      table_name = "sample_info",
      service_account_json_path = service_account_json_path
  }

  call CreateTables as CreatePetTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "pet",
      max_table_id = GetSampleIds.max_table_id,
      schema = pet_schema,
      superpartitioned = "true",
      partitioned = "true",
      uuid = "",
      service_account_json_path = service_account_json_path,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call CreateTables as CreateVetTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "vet",
      max_table_id = GetSampleIds.max_table_id,
      schema = vet_schema,
      superpartitioned = "true",
      partitioned = "true",
      uuid = "",
      service_account_json_path = service_account_json_path,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call CheckForDuplicateData {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      sample_names = external_sample_names,
      service_account_json_path = service_account_json_path,
      output_directory = output_directory,
      run_uuid = SetLock.run_uuid
  }

  call CreateFOFNs {
    input:
        input_vcf_list = write_lines(input_vcfs),
        sample_name_list = write_lines(external_sample_names),
        batch_size = batch_size,
        run_uuid = SetLock.run_uuid
  }

  scatter (i in range(length(CreateFOFNs.vcf_batch_fofns))) {
    call CreateImportTsvs {
      input:
        input_vcfs = read_lines(CreateFOFNs.vcf_batch_fofns[i]),
        sample_names = read_lines(CreateFOFNs.vcf_sample_name_fofns[i]),
        interval_list = interval_list,
        service_account_json_path = service_account_json_path,
        sample_map = select_first([GetSampleIds.sample_map, sample_map]),
        drop_state = drop_state,
        drop_state_includes_greater_than = drop_state_includes_greater_than,
        output_directory = output_directory,
        gatk_override = gatk_override,
        docker = docker_final,
        preemptible_tries = preemptible_tries,
        run_uuid = SetLock.run_uuid,
        duplicate_check_passed = CheckForDuplicateData.done
    }
  }

  scatter (i in range(GetSampleIds.max_table_id)) {
    call LoadTable as LoadPetTable {
    input:
      project_id = project_id,
      table_id = i + 1,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "pet",
      superpartitioned = "true",
      schema = pet_schema,
      service_account_json_path = service_account_json_path,
      table_creation_done = CreatePetTables.done,
      tsv_creation_done = CreateImportTsvs.done,
      docker = docker_final,
      run_uuid = SetLock.run_uuid
    }
  }

  scatter (i in range(GetSampleIds.max_table_id)) {
    call LoadTable as LoadVetTable {
    input:
      project_id = project_id,
      table_id = i + 1,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "vet",
      superpartitioned = "true",
      schema = vet_schema,
      service_account_json_path = service_account_json_path,
      table_creation_done = CreateVetTables.done,
      tsv_creation_done = CreateImportTsvs.done,
      docker = docker_final,
      run_uuid = SetLock.run_uuid
    }
  }

  call SetIsLoadedColumn {
    input:
      load_vet_done = LoadVetTable.done,
      load_pet_done = LoadPetTable.done,
      dataset_name = dataset_name,
      gvs_ids = GetSampleIds.gvs_ids,
      service_account_json_path = service_account_json_path,
      project_id = project_id,
      preemptible_tries = preemptible_tries
  }

  call ReleaseLock {
    input:
      # run_uuid = SetLock.run_uuid,
      output_directory = output_directory,
      load_sample_info_done = SetIsLoadedColumn.done,
      # load_pet_done = LoadPetTable.done,
      # load_vet_done = LoadVetTable.done,
      service_account_json_path = service_account_json_path,
      preemptible_tries = preemptible_tries
  }

  output {
    Boolean loaded_in_gvs = true
  }
}

# we create a lock file in the output directory with a uuid for this run of ImportGenomes.
# other tasks (TSV creation, bq load) check that the lock file exists and contains the run_uuid
# specific to this task.
task SetLock {
  meta {
    volatile: true
  }

  input {
    String output_directory
    String? service_account_json_path

    # runtime
    Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    # generate uuid for this run
    RUN_UUID=$(dbus-uuidgen)
    echo $RUN_UUID | tee RUN_UUID_STRING

    DIR="~{output_directory}/"

    # check for existing lock file
    LOCKFILE="LOCKFILE"
    HAS_LOCKFILE=$(gsutil ls "${DIR}${LOCKFILE}" | wc -l)
    if [ $HAS_LOCKFILE -gt 0 ]; then
      echo "ERROR: lock file in place. Check whether another run of ImportGenomes with this output directory is in progress or a previous run had an error.
            If you would like to proceed, run the following command and re-run the workflow: \
            gsutil rm ${DIR}${LOCKFILE} \
            " && exit 1
    else  # put the lock file in place
      echo "Setting lock file with UUID ${RUN_UUID}"
      echo $RUN_UUID > $LOCKFILE
      gsutil cp $LOCKFILE "${DIR}${LOCKFILE}" || { echo "Error uploading lockfile to ${DIR}${LOCKFILE}" 1>&2 ; exit 1; }
    fi
  >>>

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }

  output {
    String run_uuid = read_string("RUN_UUID_STRING")
  }
}

task ReleaseLock {
  meta {
    volatile: true
  }

  input {
    # String run_uuid
    String output_directory
    File load_sample_info_done
    # Array[String] load_pet_done
    # Array[String] load_vet_done
    String? service_account_json_path

    # runtime
    Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi


    LOCKFILE="~{output_directory}/LOCKFILE"
    EXISTING_LOCK_ID=$(gsutil cat ${LOCKFILE})

    gsutil rm $LOCKFILE

    # if [ ${EXISTING_LOCK_ID} = ${CURRENT_RUN_ID} ]; then
    #   gsutil rm $LOCKFILE
    # else
    #   echo "ERROR: found mismatched lockfile containing run ${EXISTING_LOCK_ID}, which does not match this run ${CURRENT_RUN_ID}." 1>&2
    #   exit 1
    # fi
  >>>

    runtime {
      docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
      memory: "1 GB"
      disks: "local-disk 10 HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
    }
}


task CheckForDuplicateData {
    input{
      String project_id
      String dataset_name
      Array[String] sample_names
      String? service_account_json_path
      # needed only for lockfile
      String output_directory
      String run_uuid
      # runtime
      Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'


  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # check for existence of the correct lockfile
    LOCKFILE="~{output_directory}/LOCKFILE"
    EXISTING_LOCK_ID=$(gsutil cat ${LOCKFILE}) || { echo "Error retrieving lockfile from ${LOCKFILE}" 1>&2 ; exit 1; }
    CURRENT_RUN_ID="~{run_uuid}"

    if [ "${EXISTING_LOCK_ID}" != "${CURRENT_RUN_ID}" ]; then
      echo "ERROR: found mismatched lockfile containing run ${EXISTING_LOCK_ID}, which does not match this run ${CURRENT_RUN_ID}." 1>&2
      exit 1
    fi

    INFO_SCHEMA_TABLE="~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS"
    touch duplicates
    cat ~{write_lines(sample_names)} | sort > sorted_names.txt

    # Check that the table has been created yet
    set +e
    bq show --project_id ~{project_id} $TABLE > /dev/null
    BQ_SHOW_RC=$?
    set -e
    if [ $BQ_SHOW_RC -eq 0 ]; then
      bq --location=US --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
        "SELECT s.sample_name FROM ${INFO_SCHEMA_TABLE} p JOIN ~{dataset_name}.sample_info s ON (p.partition_id = CAST(s.sample_id AS STRING)) WHERE table_name like 'pet_%'" | \
        sed -e '/sample_name/d' | sort > bq_sorted_names.txt
      comm -12 sorted_names.txt bq_sorted_names.txt > duplicates
    fi

    # true if there is data in results
    if [ -s duplicates ]; then
      echo "ERROR: Trying to load samples that have already been loaded"
      cat duplicates
      exit 1
    fi

  >>>
  runtime {
      docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
      memory: "1 GB"
      disks: "local-disk 10 HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
  output {
      File done = "sorted_names.txt"
  }
}

task CreateFOFNs {
    input {
        File input_vcf_list
        File sample_name_list
        Int batch_size
        String run_uuid
    }

     command {
         set -e

         split -d -a 5 -l ~{batch_size} ~{input_vcf_list} batched_vcfs.
         split -d -a 5 -l ~{batch_size} ~{sample_name_list} batched_sample_names.
     }

     runtime {
         docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
         bootDiskSizeGb: 15
         memory: "3 GB"
         disks: "local-disk 10 HDD"
         preemptible: 3
         cpu: 1
     }

     output {
         Array[File] vcf_batch_fofns = glob("batched_vcfs.*")
         Array[File] vcf_sample_name_fofns = glob("batched_sample_names.*")
     }
}

task CreateImportTsvs {
  input {
    Array[File] input_vcfs
    Array[String] sample_names
    File interval_list
    String output_directory
    File sample_map
    String? service_account_json_path
    String? drop_state
    Boolean? drop_state_includes_greater_than = false

    Boolean call_cache_tsvs = false
    File duplicate_check_passed

    # runtime
    Int? preemptible_tries
    File? gatk_override
    String docker

    String? for_testing_only
    String run_uuid
  }

  Int disk_size = if defined(drop_state) then 50 else 75

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  meta {
    description: "Creates a tsv file for import into BigQuery"
    volatile: true
  }

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }
  }
  command <<<
      set -e

      # workaround for https://github.com/broadinstitute/cromwell/issues/3647
      export TMPDIR=/tmp

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
      ~{for_testing_only}

      if [ ~{has_service_account_file} = 'true' ]; then
          gsutil cp ~{service_account_json_path} local.service_account.json
          export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
          gcloud auth activate-service-account --key-file=local.service_account.json
      fi


      # check for existence of the correct lockfile
      LOCKFILE="~{output_directory}/LOCKFILE"
      EXISTING_LOCK_ID=$(gsutil cat ${LOCKFILE}) || { echo "Error retrieving lockfile from ${LOCKFILE}" 1>&2 ; exit 1; }
      CURRENT_RUN_ID="~{run_uuid}"

      if [ ${EXISTING_LOCK_ID} != ${CURRENT_RUN_ID} ]; then
        echo "ERROR: found mismatched lockfile containing run ${EXISTING_LOCK_ID}, which does not match this run ${CURRENT_RUN_ID}." 1>&2
        exit 1
      fi

      # translate WDL arrays into BASH arrays
      VCFS_ARRAY=(~{sep=" " input_vcfs})
      SAMPLE_NAMES_ARRAY=(~{sep=" " sample_names})

      # loop over the BASH arrays (See https://stackoverflow.com/questions/6723426/looping-over-arrays-printing-both-index-and-value)
      for i in "${!VCFS_ARRAY[@]}"; do
          input_vcf="${VCFS_ARRAY[$i]}"
          input_vcf_basename=$(basename $input_vcf)
          updated_input_vcf=$input_vcf
          input_vcf_index="${VCFS_ARRAY[$i]}.tbi"
          sample_name="${SAMPLE_NAMES_ARRAY[$i]}"

          if [ ~{has_service_account_file} = 'true' ]; then
              gsutil cp $input_vcf .
              gsutil cp $input_vcf_index .
              updated_input_vcf=$input_vcf_basename
          fi

          # check whether these files have already been generated
          DO_TSV_GENERATION='true'
          if [ ~{call_cache_tsvs} = 'true' ]; then
            echo "Checking for files to call cache"

            declare -a TABLETYPES=("sample_info" "pet" "vet")
            ALL_FILES_EXIST='true'
            for TABLETYPE in ${TABLETYPES[@]}; do
                FILEPATH="~{output_directory}/${TABLETYPE}_tsvs/**${TABLETYPE}_*_${input_vcf_basename}.tsv"
                # output 1 if no file is found
                result=$(gsutil ls $FILEPATH || echo 1)

                if [ $result == 1 ]; then
                  echo "A file matching ${FILEPATH} does not exist"
                  ALL_FILES_EXIST='false'
                else
                  if [[ $result = *"/done/"* ]]; then
                    echo "File ${FILENAME} seems to have been processed already; found at ${result}"
                    echo "Something is very wrong!"
                    exit 1
                  elif [[ $result = "~{output_directory}/${TABLETYPE}_tsvs/set_"* ]]; then
                    FILENAME=$(basename $result)
                    echo "File ${FILENAME} is in a set directory. Moving out of set directory to ~{output_directory}/${TABLETYPE}_tsvs/${FILENAME}"
                    gsutil mv $result "~{output_directory}/${TABLETYPE}_tsvs/"
                  fi
                fi
            done

            if [ $ALL_FILES_EXIST = 'true' ]; then
                DO_TSV_GENERATION='false'
                echo "Skipping TSV generation for input file ${input_vcf_basename} because the output TSV files already exist."
            fi
          fi

          if [ $DO_TSV_GENERATION = 'true' ]; then
              echo "Generating TSVs for input file ${input_vcf_basename}"

              # TODO in future when we pass the source path or gvs_id as an arg here, use a version of that for the filename too
              gatk --java-options "-Xmx7000m" CreateVariantIngestFiles \
                -V ${updated_input_vcf} \
                -L ~{interval_list} \
                ~{"-IG " + drop_state} \
                --ignore-above-gq-threshold ~{drop_state_includes_greater_than} \
                --mode GENOMES \
                -SN $sample_name \
                -SNM ~{sample_map} \
                --ref-version 38

              gsutil -m mv sample_info_*.tsv ~{output_directory}/sample_info_tsvs/
              gsutil -m mv pet_*.tsv ~{output_directory}/pet_tsvs/
              gsutil -m mv vet_*.tsv ~{output_directory}/vet_tsvs/
          fi
      done

  >>>
  runtime {
      docker: docker
      memory: "3.75 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
  output {
      File done = "LOCKFILE"
  }
}

# Creates all the tables necessary for the LoadData operation
task CreateTables {
  meta {
    volatile: true
  }

	input {
      String project_id
      String dataset_name
      String datatype
      Int max_table_id
      String? schema
      String superpartitioned
      String partitioned
      String uuid
      String? service_account_json_path

      # runtime
      Int? preemptible_tries
      String docker
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    PREFIX=""
    if [ -n "~{uuid}" ]; then
      PREFIX="~{uuid}_"
    fi

    for TABLE_ID in $(seq 1 ~{max_table_id}); do
      PARTITION_STRING=""
      CLUSTERING_STRING=""
      if [ ~{partitioned} == "true" ]; then
        # assume clustering as well
        let "PARTITION_START=(${TABLE_ID}-1)*4000+1"
        let "PARTITION_END=$PARTITION_START+4000"
        let "PARTITION_STEP=1"
        PARTITION_FIELD="sample_id"
        CLUSTERING_FIELD="location"
        PARTITION_STRING="--range_partitioning=$PARTITION_FIELD,$PARTITION_START,$PARTITION_END,$PARTITION_STEP"
        CLUSTERING_STRING="--clustering_fields=$CLUSTERING_FIELD"
      fi

      if [ ~{superpartitioned} = "true" ]; then
        printf -v PADDED_TABLE_ID "%03d" ${TABLE_ID}
        TABLE="~{dataset_name}.${PREFIX}~{datatype}_${PADDED_TABLE_ID}"
      else
        TABLE="~{dataset_name}.${PREFIX}~{datatype}"
      fi

      # Check that the table has not been created yet
      set +e
      bq show --project_id ~{project_id} $TABLE > /dev/null
      BQ_SHOW_RC=$?
      echo $BQ_SHOW_RC > done.txt
      set -e
      if [ $BQ_SHOW_RC -ne 0 ]; then
        echo "making table $TABLE"
        bq --location=US mk ${PARTITION_STRING} ${CLUSTERING_STRING} --project_id=~{project_id} $TABLE ~{schema}
      fi
    done
  >>>

  output {
    File done = "done.txt"
  }

  runtime {
    docker: docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }
}

task LoadTable {
  meta {
    volatile: true
  }

  input {
    String project_id
    String table_id
    String dataset_name
    String storage_location
    String datatype
    String superpartitioned
    String? schema
    String? service_account_json_path
    File table_creation_done
    Array[File] tsv_creation_done
    String run_uuid

    String docker
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    DIR="~{storage_location}/~{datatype}_tsvs/"
    # check for existence of the correct lockfile
    LOCKFILE="~{storage_location}/LOCKFILE"
    EXISTING_LOCK_ID=$(gsutil cat ${LOCKFILE}) || { echo "Error retrieving lockfile from ${LOCKFILE}" 1>&2 ; exit 1; }
    CURRENT_RUN_ID="~{run_uuid}"

    if [ "${EXISTING_LOCK_ID}" != "${CURRENT_RUN_ID}" ]; then
    echo "ERROR: found mismatched lockfile containing run ${EXISTING_LOCK_ID}, which does not match this run ${CURRENT_RUN_ID}." 1>&2
    exit 1
    fi

    DIR="~{storage_location}/~{datatype}_tsvs/"

    printf -v PADDED_TABLE_ID "%03d" ~{table_id}

    FILES="~{datatype}_${PADDED_TABLE_ID}_*"

    NUM_FILES=$(gsutil ls "${DIR}${FILES}" | wc -l)

    if [ ~{superpartitioned} = "true" ]; then
      TABLE="~{dataset_name}.${PREFIX}~{datatype}_${PADDED_TABLE_ID}"
    else
      TABLE="~{dataset_name}.${PREFIX}~{datatype}"
    fi

    if [ $NUM_FILES -gt 0 ]; then
        # get list of of pet files and their byte sizes
        echo "Getting file sizes(bytes), paths to each file, and determining sets for chunking."
        echo -e "bytes\tfile_path\tsum_bytes\tset_number" > ~{datatype}_du_sets.txt
        # tr to replace each space -> tab, squeeze (remove) "empty" tabs,
        gsutil du "${DIR}${FILES}" | tr " " "\t" | tr -s "\t" | sed "/~{datatype}_tsvs\/$/d" | awk '{s+=$1}{print $1"\t"$2"\t"s"\t" (1+int(s / 16000000000000))}' >> ~{datatype}_du_sets.txt

        # per set, load table
        for set in $(sed 1d ~{datatype}_du_sets.txt | cut -f4 | sort | uniq)
        do
          # move set data into separate directory
          echo "Moving set $set data into separate directory."
          awk -v set="$set" '$4 == set' ~{datatype}_du_sets.txt | cut -f2 | gsutil -m mv -I "${DIR}set_${set}/" 2> copy.log

          # execute bq load command, get bq load job id, add details per set to tmp file
          echo "Running BigQuery load for set $set."
          bq load --quiet --nosync --location=US --project_id=~{project_id} --skip_leading_rows=1 --source_format=CSV -F "\t" \
            "$TABLE" "${DIR}set_${set}/${FILES}" ~{schema} > status_bq_submission

          cat status_bq_submission | tail -n 1 > status_bq_submission_last_line

          bq_job_id=$(sed 's/.*://' status_bq_submission_last_line)

          echo $bq_job_id
          # add job ID as key and gs path to the data set uploaded as value
          echo -e "${bq_job_id}\t${set}\t${DIR}set_${set}/" >> bq_load_details.tmp
        done

        # for each bq load job submitted, run bq wait, capture success/failure to tmp file
        while IFS="\t" read -r line_bq_load
        do
          bq wait --project_id=~{project_id} $(echo "$line_bq_load" | cut -f1) > bq_wait_status

          # capture SUCCESS or FAILURE, echo to file
          wait_status=$(sed '6q;d' bq_wait_status | tr " " "\t" | tr -s "\t" | cut -f3)
          echo "$wait_status" >> bq_wait_details.tmp
        done < bq_load_details.tmp

        # combine load status and wait status into final report
        paste bq_load_details.tmp bq_wait_details.tmp > bq_final_job_statuses.txt

        # move files from each set into set-level "done" directories
        for set in $(sed 1d ~{datatype}_du_sets.txt | cut -f4 | sort | uniq)
        do
          # move files from each set into set-level "done" directories
          echo "Moving set $set data into done directory."
          gsutil -m mv "${DIR}set_${set}/${FILES}" "${DIR}set_${set}/done/" >> gsutil_mv_done.log
        done

    else
        echo "no ${FILES} files to process in $DIR"
    fi
  >>>

  runtime {
    docker: docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 0
    cpu: 1
  }

  output {
    File done = "bq_final_job_statuses.txt"
    File? manifest_file = "~{datatype}_du_sets.txt"
    File? final_job_statuses = "bq_final_job_statuses.txt"
    File? mv_log = "gsutil_mv_done.log"
  }
}

task TerminateWorkflow {
  input {
    String message
  }

  command <<<
      set -e
      echo ~{message}
      exit 1
  >>>
  runtime {
      docker: "python:3.8-slim-buster"
      memory: "1 GB"
      disks: "local-disk 10 HDD"
      preemptible: 3
      cpu: 1
  }
}

task SetIsLoadedColumn {
  meta {
    volatile: true
  }

  input {
    Array[File] load_vet_done
    Array[File] load_pet_done
    String dataset_name
    String project_id
    File gvs_ids
    String? service_account_json_path
    Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  Array[String] gvs_id_array = read_lines(gvs_ids)

  command <<<
    set -ex

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # set is_loaded to true if there is a corresponding pet table partition with rows for that sample_id
    bq --location=US --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
    "UPDATE ~{dataset_name}.sample_info SET is_loaded = true WHERE sample_id IN (SELECT CAST(partition_id AS INT64) from ~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS WHERE partition_id in ('~{sep="\',\'" gvs_id_array}') AND total_logical_bytes > 0 AND table_name LIKE \"pet_%\")"

    # output updated schema into file for output
    bq show --project_id=~{project_id} ~{dataset_name}.sample_info > sample_info_schema.txt
>>>

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }

  output {
    File done = "sample_info_schema.txt"
  }
}

task GetSampleIds {
  meta {
    volatile: true
  }

  input {
    Array[String] external_sample_names
    String project_id
    String dataset_name
    String table_name
    String? service_account_json_path
    Int samples_per_table = 4000

    # runtime
    Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  Int num_samples = length(external_sample_names)

  command <<<
      set -ex
      if [ ~{has_service_account_file} = 'true' ]; then
        gsutil cp ~{service_account_json_path} local.service_account.json
        gcloud auth activate-service-account --key-file=local.service_account.json
      fi

      echo "project_id = ~{project_id}" > ~/.bigqueryrc

      # get the current maximum id, or 0 if there are none
      bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
        "SELECT IFNULL(MIN(sample_id),0) as min, IFNULL(MAX(sample_id),0) as max FROM ~{dataset_name}.~{table_name} where sample_name in ('~{sep="\',\'" external_sample_names}')" > results

      # prep for being able to return min table id
      min_sample_id=$(tail -1 results | cut -d, -f1)
      max_sample_id=$(tail -1 results | cut -d, -f2)

      python3 -c "from math import ceil; print(ceil($max_sample_id/~{samples_per_table}))" > max_sample_id
      python3 -c "from math import ceil; print(ceil($min_sample_id/~{samples_per_table}))" > min_sample_id

      bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false -n ~{num_samples} \
        "SELECT sample_id, sample_name FROM ~{dataset_name}.~{table_name} where sample_name in ('~{sep="\',\'" external_sample_names}')" > sample_map

      cut -d, -f1 sample_map > gvs_ids

  >>>
  runtime {
      docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
      memory: "1 GB"
      disks: "local-disk 10 HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
  output {
      Int max_table_id = ceil(read_float("max_sample_id"))
      Int min_table_id = ceil(read_float("min_sample_id"))
      File sample_map = "sample_map"
      File gvs_ids = "gvs_ids"
  }
}
