package org.broadinstitute.hellbender.tools.funcotator;


import org.apache.http.HttpEntity;
import org.apache.http.HttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.*;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.ArrayList;


/**
 * Class to copy a file using {@link CloseableHttpClient}.
 * Operates using paths.
 *
 * Created by Hailey on 8/5/21.
 */
public class FuncotatorDataSourceBundlerHttpClient {

    //==================================================================================================================
    // Standard logger:
    private final static Logger logger = LogManager.getLogger(FuncotatorDataSourceBundlerHttpClient.class);

    //==================================================================================================================
    // Public Static Members:
    public static final String ENSEMBL_CONFIG_NAME          = "ensembl.config";
    public static final String MANIFEST_FILE_NAME           = "MANIFEST.txt";
    public static final String TEMPLATE_CONFIG_FILE_NAME    = "template.config";
    public static final String README_FILE_NAME             = "README.txt";
    public static final String SCRIPT_PATH                  = "./scripts/funcotator/data_sources/fixGencodeOrdering.py";

    //==================================================================================================================
    // Private Static Members:
    private static final int BUFFER_SIZE_BYTES    = 1024 * 1024;

    //==================================================================================================================
    // Private Members:

    // Data variables:
    protected Path outputFolder;

    protected String dsOrganism;
    protected String fileName;
    protected String fastaFileName;
    protected String dsURL;
    protected Path dsPath;
    protected Path dsUnzipPath;
    protected String dsFastaURL;
    protected Path dsFastaPath;
    protected Path   dsFastaUnzipPath;
    protected Path   gtfIndexFilePath;
    protected String baseURL;
    protected String speciesName;
    protected String baseFastaURL;
    protected Path outputDestination;
    protected Path outputUnzippedDest;
    protected Path outputFastaDest;
    protected Path outputFastaUnzipDest;
    protected Path outputIndexDest;
    protected Path configFilePath;
    protected Path metadataFilePath;
    protected String dsGtfReadMeURL;
    protected String dsFastaReadMeURL;
    protected Path dsGtfReadMePath;
    protected Path dsFastaReadMePath;
    protected Path dsFastaDictPath;
    protected Path dsReorderedGtfPath;

    // Copy buffer:
    public static byte[] copyBuffer = new byte[BUFFER_SIZE_BYTES];

    //==================================================================================================================
    // Constructors:

    /**
     * {@link FuncotatorDataSourceBundlerHttpClient}
     * This internal constructor is to be used by the class itself.
     * @param outputFolder The {@link Path} into which to place the new data source supporting files.
     * @param dsOrganism The {@link String} representing the chosen organism.
     * @param speciesName The {@link String} representing the chosen division.
     * @param baseURL The {@link String} representing the base url for the chosen organism.
     * @param baseFastaURL The {@link String} representing the base url for the fasta file for the chosen organism.
     */
    protected FuncotatorDataSourceBundlerHttpClient(final Path outputFolder, final String dsOrganism, final String speciesName, final String baseURL, final String baseFastaURL) {

        this.outputFolder  = outputFolder;

        this.dsOrganism             = dsOrganism;
        this.speciesName            = speciesName;
        this.baseURL                = baseURL;
        this.baseFastaURL           = baseFastaURL;

        this.fileName               = FuncotatorDataSourceBundlerUtils.getDatasourceBaseName(this.dsOrganism, this.speciesName);
        this.fastaFileName          = FuncotatorDataSourceBundlerUtils.getFastaFileName(this.dsOrganism, this.speciesName);

        this.dsURL                  = baseURL + speciesName + "/" + fileName + "." + DataSourceUtils.GTF_GZ_EXTENSION;
        this.dsFastaURL             = baseFastaURL + speciesName + "/" + DataSourceUtils.CDNA_EXTENSION + fastaFileName + "." + DataSourceUtils.FASTA_GZ_EXTENSION;

        this.dsPath                 = IOUtils.getPath(outputFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fileName + "." + DataSourceUtils.GTF_GZ_EXTENSION);
        this.dsUnzipPath            = IOUtils.getPath(outputFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fileName + DataSourceUtils.GTF_UNZIPPED_EXTENSION);

        this.dsFastaPath            = IOUtils.getPath(outputFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fastaFileName + "." + DataSourceUtils.FASTA_GZ_EXTENSION);
        this.dsFastaUnzipPath       = IOUtils.getPath(outputFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fastaFileName + DataSourceUtils.FASTA_UNZIPPED_EXTENSION);

        this.gtfIndexFilePath       = IOUtils.getPath(outputFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fileName + "." + DataSourceUtils.GTF_GZ_EXTENSION + DataSourceUtils.IDX_EXTENSION);

        this.outputDestination      = this.dsPath.toAbsolutePath();
        this.outputUnzippedDest     = this.dsUnzipPath.toAbsolutePath();
        this.outputFastaDest        = this.dsFastaPath.toAbsolutePath();
        this.outputFastaUnzipDest   = this.dsFastaUnzipPath.toAbsolutePath();
        this.outputIndexDest        = this.gtfIndexFilePath.toAbsolutePath();

        this.configFilePath         = IOUtils.getPath(outputFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + ENSEMBL_CONFIG_NAME);
        this.metadataFilePath       = IOUtils.getPath(outputFolder + "/");

        this.dsGtfReadMeURL         = baseURL + speciesName + "/" + DataSourceUtils.README_EXTENSION;
        this.dsFastaReadMeURL       = baseURL + speciesName + "/" + DataSourceUtils.CDNA_EXTENSION + DataSourceUtils.README_EXTENSION;

        this.dsGtfReadMePath        = IOUtils.getPath(outputFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + DataSourceUtils.GTF_README_EXTENSION);
        this.dsFastaReadMePath      = IOUtils.getPath(outputFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + DataSourceUtils.FASTA_README_EXTENSION);
        this.dsFastaDictPath        = IOUtils.getPath(outputFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fileName + DataSourceUtils.FASTA_DICT_EXTENSION);
        this.dsReorderedGtfPath     = IOUtils.getPath(outputFolder + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fileName + DataSourceUtils.REORDERED_EXTENSION + DataSourceUtils.GTF_UNZIPPED_EXTENSION);
    }

    /**
     * Extract the Gzipped files that were downloaded for these datasources.
     * @param doOverwrite If {@code True} will overwrite output files.
     *                    Otherwise will throw an exception if output files already exist.
     */
    public void extractGzippedFiles(final boolean doOverwrite) {
        FuncotatorDataSourceBundlerUtils.extractGzFile(outputDestination.toString(), dsUnzipPath.toString(), doOverwrite);
        FuncotatorDataSourceBundlerUtils.extractGzFile(outputFastaDest.toString(), dsFastaUnzipPath.toString(), doOverwrite);
    }

    /**
     * Download all remote files required to create datasources for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public void downloadDataSources() {
        // Download the gtf file:
        downloadFile(dsURL, outputDestination);

        // Download the fasta file:
        downloadFile(dsFastaURL, outputFastaDest);

        // Download the gtf ReadMe file for specific data source file:
        downloadFile(dsGtfReadMeURL, dsGtfReadMePath);

        // Download the fasta ReadMe file for specific data source file:
        downloadFile(dsFastaReadMeURL, dsFastaReadMePath);
    }

    /**
     * Build a config file for the data source we have downloaded.
     */
    public void buildConfigFile() {
        try ( FileWriter writer = new FileWriter(configFilePath.toAbsolutePath().toString());
              BufferedWriter buffer = new BufferedWriter(writer) )
        {
            buffer.write(
                "name = Ensembl\n" +
                "version = 104\n" +
                "src_file = " + fileName + ".REORDERED.gtf"+ "\n" +
                "origin_location = " + dsURL + " \n" +
                "preprocessing_script = FuncotatorDataSourceBundler \n" +
                "\n" +
                "# Supported types:\n" +
                "# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID\n" +
                "# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location\n" +
                "# gencode      -- Custom datasource class for GENCODE\n" +
                "#\tcosmic       -- Custom datasource class for COSMIC\n" +
                "type = gencode\n" +
                "\n" +
                "# Required field for GENCODE files.\n" +
                "# Path to the FASTA file from which to load the sequences for GENCODE transcripts:\n" +
                "gencode_fasta_path = " + fastaFileName + ".fasta" + "\n" +
                "\n" +
                "# Required field for simpleXSV files.\n" +
                "# Valid values:\n" +
                "#     GENE_NAME\n" +
                "#     TRANSCRIPT_ID\n" +
                "xsv_key = GENE_NAME\n" +
                "\n" +
                "# Required field for simpleXSV files.\n" +
                "# The 0-based index of the column containing the key on which to match\n" +
                "xsv_key_column = 0\n" +
                "\n" +
                "# Required field for simpleXSV AND locatableXSV files.\n" +
                "# The delimiter by which to split the XSV file into columns.\n" +
                "xsv_delimiter = ,\n" +
                "\n" +
                "# Required field for simpleXSV files.\n" +
                "# Whether to permissively match the number of columns in the header and data rows\n" +
                "# Valid values:\n" +
                "#     true\n" +
                "#     false\n" +
                "xsv_permissive_cols = true\n" +
                "\n" +
                "# Required field for locatableXSV files.\n" +
                "# The 0-based index of the column containing the contig for each row\n" +
                "contig_column =\n" +
                "\n" +
                "# Required field for locatableXSV files.\n" +
                "# The 0-based index of the column containing the start position for each row\n" +
                "start_column =\n" +
                "\n" +
                "# Required field for locatableXSV files.\n" +
                "# The 0-based index of the column containing the end position for each row\n" +
                "end_column =\n"
            );
        } catch (final IOException e) {
            throw new UserException("Error. Unable to build file: " + configFilePath);
        }
    }

    /**
     * Build the template config file in the correct folder.
     */
    public void buildTemplateConfigFile() {
        try ( FileWriter writer = new FileWriter(metadataFilePath.toAbsolutePath() + "/" + TEMPLATE_CONFIG_FILE_NAME);
              BufferedWriter buffer = new BufferedWriter(writer) ) {

            buffer.write(
            "name = Achilles\n" +
                "version = 110303\n" +
                "src_file = achilles_lineage_results.import.txt\n" +
                "origin_location = UNKNOWN\n" +
                "preprocessing_script =\n" +
                "\n" +
                "# Supported types:\n" +
                "# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID\n" +
                "# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location\n" +
                "# gencode      -- Custom datasource class for GENCODE\n" +
                "# cosmic       -- Custom datasource class for COSMIC\n" +
                "type = simpleXSV\n" +
                "\n" +
                "# Required field for GENCODE files.\n" +
                "# Path to the FASTA file from which to load the sequences for GENCODE transcripts:\n" +
                "gencode_fasta_path =\n" +
                "\n" +
                "# Required field for simpleXSV files.\n" +
                "# Valid values:\n" +
                "#     GENE_NAME\n" +
                "#     TRANSCRIPT_ID\n" +
                "xsv_key = GENE_NAME\n" +
                "\n" +
                "# Required field for simpleXSV files.\n" +
                "# The 0-based index of the column containing the key on which to match\n" +
                "xsv_key_column = 0\n" +
                "\n" +
                "# Required field for simpleXSV AND locatableXSV files.\n" +
                "# The delimiter by which to split the XSV file into columns.\n" +
                "xsv_delimiter = ,\n" +
                "\n" +
                "# Required field for simpleXSV files.\n" +
                "# Whether to permissively match the number of columns in the header and data rows\n" +
                "# Valid values:\n" +
                "#     true\n" +
                "#     false\n" +
                "xsv_permissive_cols = true\n" +
                "\n" +
                "# Required field for locatableXSV files.\n" +
                "# The 0-based index of the column containing the contig for each row\n" +
                "contig_column =\n" +
                "\n" +
                "# Required field for locatableXSV files.\n" +
                "# The 0-based index of the column containing the start position for each row\n" +
                "start_column =\n" +
                "\n" +
                "# Required field for locatableXSV files.\n" +
                "# The 0-based index of the column containing the end position for each row\n" +
                "end_column ="
            );

        } catch (final IOException e) {
            throw new UserException("Error. Unable to make template config file in location: " + metadataFilePath + "/" + TEMPLATE_CONFIG_FILE_NAME);
        }
    }
    /**
     * Build a ReadMe file in the correct folder.
     */
    public void buildReadMeFile() {
        try ( final FileWriter writer = new FileWriter(metadataFilePath.toAbsolutePath() + "/" + README_FILE_NAME);
              final BufferedWriter buffer = new BufferedWriter(writer) )
        {
            buffer.write(
            "################################################################################\n" +
                "# Funcotator Data Sources Bundler Package README\n" +
                "################################################################################\n" +
                "\n" +
                "+---------------------------------------------+ \n" +
                "| Data Source Version Information             |\n" +
                "+---------------------------------------------+ \n" +
                "\n" +
                "Version:          0.0." + FuncotatorDataSourceBundlerUtils.getCurrentDateString() + "\n" +
                "Use Case:         species name\n" +
                "Source:           ./gatk -bundler.dsOrganism -species-name bundler.speciesName \n" +
                "Alternate Source: ./gatk -bundler.dsOrganism -species-name bundler.speciesName \n" +
                "\n" +
                "################################################################################\n" +
                "\n" +
                "+---------------------------------------------+ \n" +
                "| README                                      | \n" +
                "+---------------------------------------------+ \n" +
                "\n" +
                "This is a collection of data sources to be used in conjunction with Funcotator\n" +
                "to annotate data samples for a variety of species. \n" +
                "\n" +
                "This folder is a top-level Data Sources Folder for The Broad Institute's \n" +
                "Funcotator Data Source Bundler tool.  When running Funcotator, pass the path to this directory in\n" +
                "as a command-line argument:\n" +
                "\n" +
                "   ./gatk Funcotator --data-sources-path PATH/TO/THIS/FOLDER ...\n" +
                "\n" +
                "For more information on Funcotator, see the GATK development github site:\n" +
                "\n" +
                "\thttps://github.com/broadinstitute/gatk\n" +
                "\n" +
                "################################################################################\n" +
                "\n" +
                "+---------------------------------------------+ \n" +
                "| Data Sources                                |\n" +
                "+---------------------------------------------+ \n" +
                "\n" +
                "Using this Data Sources Folder will enable the following data sources:\n" +
                "--------------------\n" +
                "\n" +
                "  ensembl\n" +
                "--------------------\n" +
                "  The ENSEMBL Project produces high quality reference gene annotation and experimental validation for over 50,000 genomes. \n"
            );
        } catch (IOException e) {
            throw new UserException("Error. Unable to make ReadMe file in location: " + metadataFilePath + "/" + README_FILE_NAME);
        }
    }

    /**
     * Build a manifest file in the correct folder.
     */
    public void buildManifestFile() {
        try ( FileWriter writer = new FileWriter(metadataFilePath.toAbsolutePath() + "/" + MANIFEST_FILE_NAME);
              BufferedWriter buffer = new BufferedWriter(writer) )
        {
            buffer.write(
            "Version:          0.0." + FuncotatorDataSourceBundlerUtils.getCurrentDateString() + "\n" +
                "Use Case:         " + speciesName + "\n" +
                "Source:           ./gatk FuncotatorDataSourceBundler -" + dsOrganism + "-species-name " + speciesName + "\n" +
                "Alternate Source: ./gatk FuncotatorDataSourceBundler -" + dsOrganism + "-species-name " + speciesName + "\n"
            );

        } catch (IOException e) {
            throw new UserException("Error. Unable to make manifest file in location: " + metadataFilePath.toString() + "/" + MANIFEST_FILE_NAME);
        }
    }

    /**
     * Build an index file for the gtf data source file.
     */
    public void sortAndIndexGtfFile() {
        // Reorder gtf file:
        sortGtfFileByGenomicCoordinates(dsUnzipPath);

        // Index the GTF File
        // TODO: Fix this:
        IndexFeatureFile indexer = new IndexFeatureFile();
        indexer.indexGTF(dsReorderedGtfPath.toAbsolutePath(), gtfIndexFilePath.toAbsolutePath());
    }

    /**
     * Run the fixGencodeOrdering.py script to put gtf file in correct genomic coordinate order.
     * @param gtfFilePath The {@link Path} to the gtf file we want to reorder.
     */
    private void sortGtfFileByGenomicCoordinates(final Path gtfFilePath) {
        PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final List<String> args = new ArrayList<>();

        args.add(gtfFilePath.toString());
        args.add("--output-file");
        args.add(dsReorderedGtfPath.toAbsolutePath().toString());
        boolean success = executor.executeScript("./scripts/funcotator/data_sources/fixGencodeOrdering.py", null, args);
        if (!success) {
            throw new UserException("Error. Unable to sort gtf file by genomic coordinates.");
        }
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Download a file at the given URL to the given destination.
     * @param url {@link String} containing the url at which the source file is located.
     * @param dest {@link Path} representing the destination of the downloaded file.
     */
    private static void downloadFile(final String url, final Path dest) {

        logger.info("Downloading file: " + url + " -> " + dest.toUri());

        // Creating CloseableHttpClient object to access the webpage and retrieve the file:
        final CloseableHttpClient client = HttpClientBuilder.create().build();

        // Creating an HttpGet object to send the request to the server:
        final HttpGet request = new HttpGet(url);

        try {
            // Using an HttpResponse class object to catch the response from the server
            final HttpResponse response = client.execute(request);

            // The data sent by the server is obtained in this getEntity() function:
            final HttpEntity entity = response.getEntity();

            // Extracting the data from the entity object:
            try( final InputStream inputStream = entity.getContent();
                 final OutputStream outputStream = Files.newOutputStream(dest) )
            {

                // Perform the copy:
                while (true) {

                    // Read from our input:
                    final int bytesRead = inputStream.read(copyBuffer);
                    if (bytesRead == -1) {
                        break;
                    }

                    // Write to our output:
                    outputStream.write(copyBuffer, 0, bytesRead);
                }
            }
            catch (final IOException ex) {
                throw new UserException("Could not copy file: " + url + " -> " + dest.toUri(), ex);
            }
        }
        catch (final IOException ex) {
            throw new UserException("Could not obtain data from "+ url, ex);
        }
    }

    //==================================================================================================================
    // Getters / Setters:

    /**
     * @return A copy of the {@link String} used as the data source URL for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public String getDSUrl() {
        return this.dsURL;
    }

    /**
     * @return A copy of the {@link Path} used as the output destination path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getOutputDestination() {
        return this.outputDestination;
    }

    /**
     * @return A copy of the {@link Path} used as the output fasta desination path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getFastaOutputDestination() {
        return this.outputFastaDest;
    }

    /**
     * @return A copy of the {@link Path} used as the data source path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getDSPath() {
        return this.dsPath;
    }

    /**
     * @return A copy of the {@link Path} used as the unzipped data source path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getDSUnzipPath() {
        return this.dsUnzipPath;
    }

    /**
     * @return A copy of the {@link Path} used as the unzipped fasta data source path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getFastaUnzipPath() {
        return this.dsFastaUnzipPath;
    }

    /**
     * @return A copy of the {@link Path} used as the index file path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getIndexPath() {
        return this.gtfIndexFilePath;
    }
}