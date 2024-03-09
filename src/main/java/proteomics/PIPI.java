/*
 * Copyright 2016-2019 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package proteomics;

import ProteomicsLibrary.*;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.FM.FMIndex;
import proteomics.PTM.InferPTM;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.DatasetReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;


import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.sql.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class PIPI {
    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    public static final String versionStr = "2.0";
    static final boolean useXcorr = false;
    static final int minTagLenToReduceProtDb = 5;
    public static final boolean isPtmSimuTest = false;
    static int minTagLenToExtract;
    static int maxTagLenToExtract;
    static boolean nTermSpecific;
    public static double MIN_PEAK_SUM_INFER_AA;
    static double proteinCovThres;
    static int  maxNumVarPtmConsidered;

    public static void main(String[] args) throws InterruptedException, IOException {
        long startTime = System.nanoTime();
        // Process inputs
        if (args.length != 3) {
            help();
        }
        // Set parameters
        String parameterPath = args[0].trim();
        String spectraPath = args[1].trim();
        String outputDir = args[2].trim();

        logger.info("Running PIPI version {}.", versionStr);

        String dbName = null;
        String hostName = "unknown-host";
        try {
            hostName = InetAddress.getLocalHost().getHostName();
            logger.info("Computer: {}.", hostName);
        } catch (UnknownHostException ex) {
            logger.warn("Cannot get the computer's name.");
        }

        try {
            logger.info("Spectra: {}, parameter: {}.", spectraPath, parameterPath);

            dbName = String.format(Locale.US, "PIPI.temp.db");
            new PIPI(parameterPath, spectraPath, dbName, hostName, outputDir);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
        } finally {
            if (dbName != null) {
                (new File(dbName)).delete();
                (new File(dbName + "-wal")).delete();
                (new File(dbName + "-shm")).delete();
            }
        }

        double totalHour = (double) (System.nanoTime() - startTime) * 1e-9 / 3600;
        logger.info("Running time: {} hours.", totalHour);
        logger.info("Done!");
    }

    private PIPI(String parameterPath, String spectraPath, String dbName, String hostName, String outputDir) throws Exception {
        // Get the parameter map
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double ms1Tolerance = Double.valueOf(parameterMap.get("ms1_tolerance"));
        double leftInverseMs1Tolerance = 1 / (1 + ms1Tolerance * 1e-6);
        double rightInverseMs1Tolerance = 1 / (1 - ms1Tolerance * 1e-6);
        int ms1ToleranceUnit = Integer.valueOf(parameterMap.get("ms1_tolerance_unit"));
        double minClear = Double.valueOf(parameterMap.get("min_clear_mz"));
        double maxClear = Double.valueOf(parameterMap.get("max_clear_mz"));

        minTagLenToExtract = Integer.valueOf(parameterMap.get("minTagLenToExtract"));
        maxTagLenToExtract = Integer.valueOf(parameterMap.get("maxTagLenToExtract"));
        nTermSpecific = Boolean.valueOf(parameterMap.get("nTermSpecific"));
        MIN_PEAK_SUM_INFER_AA = Double.valueOf(parameterMap.get("MIN_PEAK_SUM_INFER_AA"));
        proteinCovThres = Double.valueOf(parameterMap.get("proteinCovThres"));
        maxNumVarPtmConsidered = Integer.valueOf(parameterMap.get("maxNumVarPtmConsidered"));

        String[] tempArray = parameterMap.get("ms_level").split(",");
        Set<Integer> msLevelSet = new HashSet<>(tempArray.length + 1, 1);
        for (String temp : tempArray) {
            msLevelSet.add(Integer.valueOf(temp));
        }

        logger.info("Loading parameters and build fmIndex...");
        BuildIndex buildIndex = new BuildIndex(parameterMap);
        MassTool massTool = buildIndex.returnMassTool();
        InferPTM inferPTM = buildIndex.getInferPTM();

        logger.info("Reading spectra...");
        File spectraFile = new File(spectraPath);
        DatasetReader datasetReader;
        JMzReader[] spectraParserArray;
        String sqlPath = "jdbc:sqlite:" + dbName;
        Class.forName("org.sqlite.JDBC").newInstance();
        Map<Integer, String> fileIdNameMap = new HashMap<>();
        Map<String, Integer> fileNameIdMap = new HashMap<>();
        if ((!spectraFile.exists())) {
            throw new FileNotFoundException("The spectra file not found.");
        }

        if ( ! spectraFile.isDirectory()) {
            spectraParserArray = new JMzReader[1];
            JMzReader spectraParser;
            String ext = spectraPath.substring(spectraPath.lastIndexOf(".")+1);
            if (ext.contentEquals("mzXML")) {
                spectraParser = new MzXMLFile(spectraFile);
            } else if (ext.toLowerCase().contentEquals("mgf")) {
                spectraParser = new MgfFile(spectraFile);
            } else {
                throw new Exception(String.format(Locale.US, "Unsupported file format %s. Currently, PIPI only support mzXML and MGF.", ext));
            }
            spectraParserArray[0] = spectraParser;
            int lastSlash = -1;
            if (spectraPath.contains("\\")) {
                lastSlash = spectraPath.lastIndexOf("\\");
            } else if (spectraPath.contains("/")){
                lastSlash = spectraPath.lastIndexOf("/");
            }
            fileIdNameMap.put(0, spectraPath.substring(lastSlash+1).split("\\.")[0].replaceAll("\\.","_"));
            fileNameIdMap.put(spectraPath.substring(lastSlash+1).split("\\.")[0].replaceAll("\\.","_"), 0);
            datasetReader = new DatasetReader(spectraParserArray, ms1Tolerance, ms1ToleranceUnit, massTool, ext, msLevelSet, sqlPath, fileIdNameMap);
        } else {
            String[] fileList = spectraFile.list(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    return name.endsWith(".mgf");
                }
            });
            spectraParserArray = new JMzReader[fileList.length];
            for (int i = 0; i < fileList.length; i++){
                spectraParserArray[i] = new MgfFile(new File(spectraPath + fileList[i]));
                fileIdNameMap.put(i, fileList[i].split("\\.")[0].replaceAll("\\.","_"));
                fileNameIdMap.put(fileList[i].split("\\.")[0].replaceAll("\\.","_"), i);

            }

            String ext = fileList[0].substring(fileList[0].lastIndexOf(".")+1);
            datasetReader = new DatasetReader(spectraParserArray, ms1Tolerance, ms1ToleranceUnit, massTool, ext, msLevelSet, sqlPath, fileIdNameMap);
        }

        SpecProcessor specProcessor = new SpecProcessor(massTool);
        Map<String, Integer> precursorChargeMap = new HashMap<>();
        Map<String, Double> precursorMassMap = new HashMap<>();
        TreeMap<Double, Set<String>> pcMassScanNameMap = new TreeMap<>();

        logger.info("Get long tags to reduce proteins...");
        int threadNum_0 = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum_0 == 0) {
            threadNum_0 = 3 + Runtime.getRuntime().availableProcessors();
        }

        Map<String, String> mgfTitleMap = new HashMap<>();
        Map<String, Integer> isotopeCorrectionNumMap = new HashMap<>();
        Map<String, Double> ms1PearsonCorrelationCoefficientMap = new HashMap<>();

        //////////==================================================
        ExecutorService threadPoolGetLongTag = Executors.newFixedThreadPool(threadNum_0);
        ArrayList<Future<GetLongTag.Entry>> taskListGetLongTag = new ArrayList<>(datasetReader.getUsefulSpectraNum() + 10);
        Connection sqlConSpecCoderX = DriverManager.getConnection(sqlPath);
        Statement sqlStateGetLongTag = sqlConSpecCoderX.createStatement();
        ResultSet sqlResSetGetLongTag = sqlStateGetLongTag.executeQuery("SELECT scanName, scanNum, precursorCharge" +
                ", precursorMass, precursorScanNo, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient FROM spectraTable");

        ReentrantLock lockGetLongTag = new ReentrantLock();

        int submitNumSpecCoderX = 0;
        Set<String> validScanSet = new HashSet<>();

        while (sqlResSetGetLongTag.next()) {
            String scanName = sqlResSetGetLongTag.getString("scanName");
            int scanNum = sqlResSetGetLongTag.getInt("scanNum");
            int precursorCharge = sqlResSetGetLongTag.getInt("precursorCharge");
            double precursorMass = sqlResSetGetLongTag.getDouble("precursorMass");
            mgfTitleMap.put(scanName, sqlResSetGetLongTag.getString("mgfTitle"));
            isotopeCorrectionNumMap.put(scanName, sqlResSetGetLongTag.getInt("isotopeCorrectionNum"));
            ms1PearsonCorrelationCoefficientMap.put(scanName, sqlResSetGetLongTag.getDouble("ms1PearsonCorrelationCoefficient"));

            int fileId = fileNameIdMap.get( scanName.split("\\.")[0] );

            precursorChargeMap.put(scanName, precursorCharge);
            precursorMassMap.put(scanName, precursorMass);
            submitNumSpecCoderX++;
            validScanSet.add(scanName);
            taskListGetLongTag.add(threadPoolGetLongTag.submit(new GetLongTag(scanNum, buildIndex, massTool, spectraParserArray[fileId], minClear, maxClear, lockGetLongTag, scanName, precursorCharge
                    , precursorMass, specProcessor, "dummyseq", ms2Tolerance)));
        }
        sqlResSetGetLongTag.close();
        sqlStateGetLongTag.close();

        int lastProgressGetLongTag = 0;
        int totalCountGetLongTag = taskListGetLongTag.size();
        int countGetLongTag = 0;
        Map<String, List<GetLongTag.TagRes>> prot_TagResList_Map = new HashMap<>();
        while (countGetLongTag < totalCountGetLongTag) {
            // record search results and delete finished ones.
            List<Future<GetLongTag.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountGetLongTag - countGetLongTag);
            for (Future<GetLongTag.Entry> task : taskListGetLongTag) {
                if (task.isDone()) {
                    if (task.get() != null ) {
                        GetLongTag.Entry entry = task.get();
                        for (String prot : entry.prot_TagResList_Map.keySet()){
                            List<GetLongTag.TagRes> tagResList = entry.prot_TagResList_Map.get(prot);
                            if  (prot_TagResList_Map.containsKey(prot)) {
                                prot_TagResList_Map.get(prot).addAll(tagResList);
                            } else {
                                List<GetLongTag.TagRes> tmpTagResList = new LinkedList<>();
                                tmpTagResList.addAll(tagResList);
                                prot_TagResList_Map.put(prot, tagResList);
                            }
                        }
                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countGetLongTag += toBeDeleteTaskList.size();
            taskListGetLongTag.removeAll(toBeDeleteTaskList);
            taskListGetLongTag.trimToSize();

            int progress = countGetLongTag * 20 / totalCountGetLongTag;
            if (progress != lastProgressGetLongTag) {
                logger.info("Getting long tags for prot {}%...", progress * 5);
                lastProgressGetLongTag = progress;
            }

            if (countGetLongTag == totalCountGetLongTag) {
                break;
            }
            Thread.sleep(6000);
        }
        // shutdown threads.
        threadPoolGetLongTag.shutdown();
        if (!threadPoolGetLongTag.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolGetLongTag.shutdownNow();
            if (!threadPoolGetLongTag.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }
        if (lockGetLongTag.isLocked()) {
            lockGetLongTag.unlock();
        }
        //   ==========END=============Get Long Tags
        //   Score for Proteins
        Map<String, Integer> protLengthMap = buildIndex.protLengthMap;
        Map<String, Double> protCoverageMap = new HashMap<>();
        for (String prot : prot_TagResList_Map.keySet()){
            int protLen = buildIndex.protSeqMap.get(prot).length();
            List<GetLongTag.TagRes> tagResList = prot_TagResList_Map.get(prot);
            Map<Integer, Double> posProductMap = new HashMap<>(protLen);
            for (int pos = 0; pos < protLen; pos++) {
                posProductMap.put(pos, 1.0); //intial, all pos is with 1
            }
            for (GetLongTag.TagRes tagRes : tagResList) {
                int relPos = tagRes.relPos;
                List<Double> normedIaaList = tagRes.normedIaaList;
                for (int aaPos = relPos; aaPos < relPos+tagRes.tagSeq.length(); aaPos++) {
                    posProductMap.put(aaPos, posProductMap.get(aaPos)*(1- normedIaaList.get(aaPos-relPos)));
                }
            } // finish aa cal
            double tempCoverage = 0;
            for (int pos : posProductMap.keySet()) {
                tempCoverage += 1-posProductMap.get(pos);
            }
            protCoverageMap.put(prot, tempCoverage/protLen);
        }

        List<Pair<String, Double>> protScoreLongList = new ArrayList<>();
        for (String prot : protCoverageMap.keySet()){
            protScoreLongList.add(new Pair<>(prot, protCoverageMap.get(prot)));
        }
        Collections.sort(protScoreLongList, Comparator.comparing(o -> o.getSecond(), Comparator.reverseOrder()));
        Set<String> reducedProtIdSet = new HashSet<>();
        int ii = 0;
        for (Pair<String, Double> pair : protScoreLongList){
            ii++;
            if (pair.getSecond() < proteinCovThres) break;
            reducedProtIdSet.add(pair.getFirst());
        }

        if (submitNumSpecCoderX < 500)
        {
            reducedProtIdSet = protLengthMap.keySet();
        }

        Iterator<String> iter = buildIndex.protSeqMap.keySet().iterator();
        while (iter.hasNext()) {
            if (!reducedProtIdSet.contains(iter.next())) {
                iter.remove();
            }
        }
        //   Get Tags and Get Peptide Candidates
        logger.info("Building decoy...");
        int threadNum_1 = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum_1 == 0) {
            threadNum_1 = 3 + Runtime.getRuntime().availableProcessors();
        }

        ExecutorService threadPoolBuildDecoyProts = Executors.newFixedThreadPool(threadNum_1);
        ArrayList<Future<BuildDecoyProts.Entry>> taskListBuildDecoyProts = new ArrayList<>(reducedProtIdSet.size() + 10);
        ReentrantLock lockSpecCoder = new ReentrantLock();

        for (String protId : reducedProtIdSet) {
            taskListBuildDecoyProts.add(threadPoolBuildDecoyProts.submit(new BuildDecoyProts(parameterMap, buildIndex, protId)));
        }

        int lastProgressBuildDecoyProts = 0;
        int totalCountBuildDecoyProts = taskListBuildDecoyProts.size();
        int countBuildDecoyProts = 0;
        double minPeptideMass = 9999;
        double maxPeptideMass = 0;

        Map<String, Set<Pair<String, Integer>>> tagProtPosMap = new HashMap<>();
        while (countBuildDecoyProts < totalCountBuildDecoyProts) {
            // record search results and delete finished ones.
            List<Future<BuildDecoyProts.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountBuildDecoyProts - countBuildDecoyProts);
            for (Future<BuildDecoyProts.Entry> task : taskListBuildDecoyProts) {
                if (task.isDone()) {
                    if (task.get() != null ) {
                        BuildDecoyProts.Entry entry = task.get();
                        String protId = entry.protId;
                        String decoyProtSeq = entry.decoyProtSeq;
                        String decoyProtId = "DECOY_"+protId;
                        buildIndex.protSeqMap.put(decoyProtId, decoyProtSeq);

                    }
                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countBuildDecoyProts += toBeDeleteTaskList.size();
            taskListBuildDecoyProts.removeAll(toBeDeleteTaskList);
            taskListBuildDecoyProts.trimToSize();

            int progress = countBuildDecoyProts * 20 / totalCountBuildDecoyProts;
            if (progress != lastProgressBuildDecoyProts) {
                logger.info("Build Decoy Prots {}%...", progress * 5);
                lastProgressBuildDecoyProts = progress;
            }

            if (countBuildDecoyProts == totalCountBuildDecoyProts) {
                break;
            }
            Thread.sleep(6000);
        }
        // shutdown threads.
        threadPoolBuildDecoyProts.shutdown();
        if (!threadPoolBuildDecoyProts.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolBuildDecoyProts.shutdownNow();
            if (!threadPoolBuildDecoyProts.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }
        if (lockSpecCoder.isLocked()) {
            lockSpecCoder.unlock();
        }

        logger.info("Generating reduced FM index...");

        BufferedWriter writerProt = new BufferedWriter(new FileWriter("catProtReduced.txt"));
        int dotPos = 0;
        int dotNum = 0;
        buildIndex.dotPosArrReduced = new int[buildIndex.protSeqMap.keySet().size()];

        for (String protId : buildIndex.protSeqMap.keySet()) {
            buildIndex.dotPosArrReduced[dotNum] = dotPos;
            buildIndex.posProtMapReduced.put(dotNum, protId);
            String protSeq = buildIndex.protSeqMap.get(protId).replace('I', 'L');
            buildIndex.protSeqMap.put(protId, protSeq);
            writerProt.write("." + protSeq.replace('I', 'L'));
            dotNum++;
            dotPos += protSeq.length()+1;
            int numOfTags = buildIndex.inferSegment.getLongTagNumForProt(protSeq);
            protLengthMap.put(protId, numOfTags);
        }
        writerProt.close();
        char[] text = buildIndex.loadFile("catProtReduced.txt", true);
        buildIndex.fmIndexReduced = new FMIndex(text);
        buildIndex.minPeptideMass = minPeptideMass;
        buildIndex.maxPeptideMass = maxPeptideMass;
        buildIndex.tagProtPosMap = tagProtPosMap;
        // writer concatenated fasta
        Map<String, String> proteinAnnotationMap;
        String dbPath = parameterMap.get("db");
        DbTool dbTool = buildIndex.dbTool;
        DbTool contaminantsDb = null;
        if (true) {
            contaminantsDb = new DbTool(null, "contaminants");
            proteinAnnotationMap = contaminantsDb.getProteinAnnotateMap();
            proteinAnnotationMap.putAll(dbTool.getProteinAnnotateMap()); // using the target annotation to replace contaminant sequence if there is conflict.
        } else {
            proteinAnnotationMap = dbTool.getProteinAnnotateMap();
        }
        BufferedWriter writer = new BufferedWriter(new FileWriter(dbPath + ".TD.fasta"));
        for (String protId : buildIndex.protSeqMap.keySet()) {
            writer.write(String.format(Locale.US, ">%s %s\n", protId, proteinAnnotationMap.getOrDefault(protId, "")));
            writer.write(buildIndex.protSeqMap.get(protId) + "\n");
        }
        writer.close();

        int threadNum = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum == 0) {
            threadNum = 3 + Runtime.getRuntime().availableProcessors();
        }
        ExecutorService threadPoolBone = Executors.newFixedThreadPool(threadNum);
        ArrayList<Future<PreSearch.Entry>> taskListBone = new ArrayList<>(validScanSet.size() + 10);
        ReentrantLock lockBone = new ReentrantLock();
        for (String scanName : validScanSet) {
            String[] scanNameStr = scanName.split("\\.");
            int precursorCharge = precursorChargeMap.get(scanName);

            double precursorMass = precursorMassMap.get(scanName);
            int scanNum = Integer.valueOf(scanNameStr[2]);
            int fileId = fileNameIdMap.get( scanNameStr[0] );
            taskListBone.add(threadPoolBone.submit(new PreSearch(scanNum, buildIndex, massTool, ms2Tolerance, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance
                    , ms1ToleranceUnit, inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass(), Math.min(precursorCharge > 1 ? precursorCharge - 1 : 1, 3)
                    , spectraParserArray[fileId], minClear, maxClear, lockBone, scanName, precursorCharge, precursorMass, specProcessor , "dasdsads", Integer.valueOf(parameterMap.get("min_peptide_length")), Integer.valueOf(parameterMap.get("max_peptide_length")))));
        }

        Map<String, PeptideInfo> allPeptideInfoMap = new HashMap<>();
        int lastProgressBone = 0;
        int totalCountBone = taskListBone.size();
        int countBone = 0;

        Map<Integer, TreeMap<Double, Set<String>>> fileId_pcMassScanNameMap = new HashMap<>();
        for (int fileId : fileIdNameMap.keySet()) {
            fileId_pcMassScanNameMap.put(fileId, new TreeMap<>());
        }
        Map<String, Peptide> scanName_TopPeptide_Map = new HashMap<>();
        Map<String, String> scanName_PepString_Map = new HashMap<>();
        Connection sqlConSpecCoder = DriverManager.getConnection(sqlPath);
        PreparedStatement sqlPreparedStatement = sqlConSpecCoder.prepareStatement("REPLACE INTO spectraTable (scanNum, scanName,  precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass, isDecoy, score, deltaLCn, deltaCn, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, otherPtmPatterns, aScore, candidates, peptideSet, whereIsTopCand, shouldPtm, hasPTM, ptmNum, isSettled) VALUES (?, ?, ?, ?,?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConSpecCoder.setAutoCommit(false);
        Map<VarPtm, Integer> varPtmCountMap = new HashMap<>();
        while (countBone < totalCountBone) {
            // record search results and delete finished ones.
            List<Future<PreSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountBone - countBone);
            for (Future<PreSearch.Entry> task : taskListBone) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        PreSearch.Entry entry = task.get();

                        for (VarPtm varPtm : entry.varPtmList){
                            if (varPtmCountMap.containsKey(varPtm)) {
                                varPtmCountMap.put(varPtm, varPtmCountMap.get(varPtm)+1);
                            }else{
                                varPtmCountMap.put(varPtm, 1);
                            }
                        }
                        allPeptideInfoMap.putAll(entry.peptideInfoMapForRef);
                        sqlPreparedStatement.setInt(1, entry.scanNum);
                        sqlPreparedStatement.setString(2, entry.scanName);
                        sqlPreparedStatement.setInt(3, entry.precursorCharge);
                        sqlPreparedStatement.setDouble(4, entry.precursorMass);
                        sqlPreparedStatement.setString(5, mgfTitleMap.get(entry.scanName));
                        sqlPreparedStatement.setInt(6, isotopeCorrectionNumMap.get(entry.scanName));
                        sqlPreparedStatement.setDouble(7, ms1PearsonCorrelationCoefficientMap.get(entry.scanName));
                        sqlPreparedStatement.setString(8, entry.labelling);
                        sqlPreparedStatement.setString(9, entry.peptide);
                        sqlPreparedStatement.setDouble(10, entry.theoMass);
                        sqlPreparedStatement.setInt(11, entry.isDecoy);
                        sqlPreparedStatement.setDouble(12, entry.score);
                        sqlPreparedStatement.setDouble(13, entry.deltaLCn);
                        sqlPreparedStatement.setDouble(14, entry.deltaCn);
                        sqlPreparedStatement.setInt(15, entry.matchedPeakNum);
                        sqlPreparedStatement.setDouble(16, entry.ionFrac);
                        sqlPreparedStatement.setDouble(17, entry.matchedHighestIntensityFrac);
                        sqlPreparedStatement.setDouble(18, entry.explainedAaFrac);
                        sqlPreparedStatement.setString(19, entry.otherPtmPatterns);
                        sqlPreparedStatement.setString(20, entry.aScore);
                        sqlPreparedStatement.setString(21, entry.candidates);
                        sqlPreparedStatement.setString(22, entry.peptideSet);
                        sqlPreparedStatement.setInt(23, -1);
                        sqlPreparedStatement.setInt(24, entry.shouldPtm);
                        sqlPreparedStatement.setInt(25, entry.hasPTM);
                        sqlPreparedStatement.setInt(26, entry.ptmNum);
                        sqlPreparedStatement.setInt(27, entry.isSettled);

                        sqlPreparedStatement.executeUpdate();
                        int fileId = fileNameIdMap.get( entry.scanName.split("\\.")[0] );
                        TreeMap<Double,Set<String>> local_pcMassScanNameMap = fileId_pcMassScanNameMap.get(fileId);
                        if (local_pcMassScanNameMap.containsKey(entry.precursorMass)) {
                            local_pcMassScanNameMap.get(entry.precursorMass).add(entry.scanName);
                        }else {
                            Set<String> scanNumSet = new HashSet<>();
                            scanNumSet.add(entry.scanName);
                            local_pcMassScanNameMap.put(entry.precursorMass, scanNumSet);
                        }
                        scanName_TopPeptide_Map.put(entry.scanName, entry.topPeptide);
                        scanName_PepString_Map.put(entry.scanName, entry.peptideSet);

                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countBone += toBeDeleteTaskList.size();
            taskListBone.removeAll(toBeDeleteTaskList);
            taskListBone.trimToSize();

            sqlConSpecCoder.commit();
            int progress = countBone * 20 / totalCountBone;
            if (progress != lastProgressBone) {
                logger.info("Pre Searching {}%...", progress * 5);
                lastProgressBone = progress;
            }

            if (countBone == totalCountBone) {
                break;
            }
            Thread.sleep(6000);
        }
        // shutdown threads.
        threadPoolBone.shutdown();
        if (!threadPoolBone.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolBone.shutdownNow();
            if (!threadPoolBone.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }
        sqlConSpecCoder.commit();
        sqlConSpecCoder.setAutoCommit(true);
        sqlConSpecCoder.close();
        if (lockBone.isLocked()) {
            lockBone.unlock();
        }

        if (threadNum == 0) {
            threadNum = 3 + Runtime.getRuntime().availableProcessors();
        }
        ExecutorService threadPoolPtm = Executors.newFixedThreadPool(threadNum);
        ArrayList<Future<PostSearch.Entry>> taskListPTM = new ArrayList<>(validScanSet.size() + 10);
        ReentrantLock lockPtm = new ReentrantLock();
        Binomial binomial = new Binomial(Integer.valueOf(parameterMap.get("max_peptide_length")) * 2);

        int n_complementaryPeptides = 5;
        for (String thisScanName : validScanSet) {
            int thisFileId = fileNameIdMap.get( thisScanName.split("\\.")[0] );
            double mass = precursorMassMap.get(thisScanName);
            int thisScanNum = Integer.valueOf( thisScanName.split("\\.")[2] );

            String thisPtmSeq = "XFAKEPEP";
            if (scanName_TopPeptide_Map.containsKey(thisScanName)) {
                thisPtmSeq = scanName_TopPeptide_Map.get(thisScanName).getVarPtmContainingSeqNow();
            }

            Map<String, TreeMap<Integer, VarPtm>> ptmSeq_posVarPtmMap_Map = new HashMap<>();
            Map<String, PeptideInfo> local_PepSeq_PeptideInfo_Map = new HashMap<>();

            TreeMap<Double, Set<String>> local_other_pcMassScanNameMap = fileId_pcMassScanNameMap.get(thisFileId);
            for (int i = 0; i <= 0; i++) {
                for (Set<String> otherScanNameSet : local_other_pcMassScanNameMap.subMap(mass + i*MassTool.PROTON - mass*ms1Tolerance*1e-6, true, mass + i*MassTool.PROTON + mass*ms1Tolerance*1e-6, true).values()) {
                    for (String otherScanName : otherScanNameSet) {
                        int otherScanNum = Integer.valueOf( otherScanName.split("\\.")[2] );
                        if (otherScanNum < thisScanNum+2000 && otherScanNum > thisScanNum-2000 && otherScanNum != thisScanNum) {
                            String otherPtmSeq = scanName_TopPeptide_Map.get(otherScanName).getVarPtmContainingSeqNow();
                            String otherFreeSeq = scanName_TopPeptide_Map.get(otherScanName).getFreeSeq();

                            if (otherPtmSeq.contentEquals(thisPtmSeq)) continue;

                            if (!ptmSeq_posVarPtmMap_Map.containsKey(otherPtmSeq)) {
                                ptmSeq_posVarPtmMap_Map.put(otherPtmSeq, scanName_TopPeptide_Map.get(otherScanName).posVarPtmResMap);
                                local_PepSeq_PeptideInfo_Map.put(otherFreeSeq, allPeptideInfoMap.get(otherFreeSeq).clone());
                            }
                        }
                    }
                }
            }

            if (ptmSeq_posVarPtmMap_Map.isEmpty()) {
                continue;
            }
            int precursorCharge = precursorChargeMap.get(thisScanName);
            double precursorMass = precursorMassMap.get(thisScanName);
            taskListPTM.add(threadPoolPtm.submit(new PostSearch(thisScanNum, buildIndex, massTool, ms1Tolerance
                    , leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, ms2Tolerance
                    , inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass(), Math.min(precursorCharge > 1 ? precursorCharge - 1 : 1, 3)
                    , spectraParserArray[thisFileId], minClear, maxClear, lockPtm, thisScanName, precursorCharge
                    , precursorMass, inferPTM, specProcessor, binomial, ptmSeq_posVarPtmMap_Map, local_PepSeq_PeptideInfo_Map)));
        }
        // check progress every minute, record results,and delete finished tasks.
        int lastProgressPtm = 0;
        int totalCountPtm = taskListPTM.size();
        int countPtm = 0;
        Connection sqlConPTM = DriverManager.getConnection(sqlPath);
        PreparedStatement sqlPreparedStatementPTM = sqlConPTM.prepareStatement("REPLACE INTO spectraTable (scanNum, scanName,  precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass, isDecoy, score, deltaLCn, deltaCn, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, otherPtmPatterns, aScore, candidates, peptideSet, whereIsTopCand, shouldPtm, hasPTM, ptmNum, isSettled) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConPTM.setAutoCommit(false);
        while (countPtm < totalCountPtm) {
            List<Future<PostSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountPtm - countPtm);
            for (Future<PostSearch.Entry> task : taskListPTM) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        PostSearch.Entry entry = task.get();

                        if (  scanName_TopPeptide_Map.containsKey(entry.scanName)
                            && entry.score < scanName_TopPeptide_Map.get(entry.scanName).getScore()) {
                            toBeDeleteTaskList.add(task);
                            continue; // if the complementary is not better dont change in sql
                        }

                        sqlPreparedStatementPTM.setInt(1, entry.scanNum);
                        sqlPreparedStatementPTM.setString(2, entry.scanName);

                        sqlPreparedStatementPTM.setInt(3, entry.precursorCharge);
                        sqlPreparedStatementPTM.setDouble(4, entry.precursorMass);
                        sqlPreparedStatementPTM.setString(5, mgfTitleMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setInt(6, isotopeCorrectionNumMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setDouble(7, ms1PearsonCorrelationCoefficientMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setString(8, entry.labelling);
                        sqlPreparedStatementPTM.setString(9, entry.peptide);
                        sqlPreparedStatementPTM.setDouble(10, entry.theoMass);
                        sqlPreparedStatementPTM.setInt(11, entry.isDecoy);
                        sqlPreparedStatementPTM.setDouble(12, entry.score);
                        sqlPreparedStatementPTM.setDouble(13, entry.deltaLCn);
                        sqlPreparedStatementPTM.setDouble(14, entry.deltaCn);
                        sqlPreparedStatementPTM.setInt(15, entry.matchedPeakNum);
                        sqlPreparedStatementPTM.setDouble(16, entry.ionFrac);
                        sqlPreparedStatementPTM.setDouble(17, entry.matchedHighestIntensityFrac);
                        sqlPreparedStatementPTM.setDouble(18, entry.explainedAaFrac);
                        sqlPreparedStatementPTM.setString(19, entry.otherPtmPatterns);
                        sqlPreparedStatementPTM.setString(20, entry.aScore);
                        sqlPreparedStatementPTM.setString(21, entry.candidates);
                        String peptideSetStr = null;
                        if (  scanName_PepString_Map.containsKey(entry.scanName) ){
                            peptideSetStr = entry.peptideSet + "," + scanName_PepString_Map.get(entry.scanName);
                        } else {
                            peptideSetStr = entry.peptideSet;
                        }
                        sqlPreparedStatementPTM.setString(22, peptideSetStr); //if a better complementary pep comes, need to insert it in the front
                        sqlPreparedStatementPTM.setInt(23, entry.whereIsTopCand);
                        sqlPreparedStatementPTM.setInt(24, entry.shouldPtm);
                        sqlPreparedStatementPTM.setInt(25, entry.hasPTM);
                        sqlPreparedStatementPTM.setInt(26, entry.ptmNum);
                        sqlPreparedStatementPTM.setInt(27, entry.isSettled);

                        sqlPreparedStatementPTM.executeUpdate();
                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countPtm += toBeDeleteTaskList.size();
            taskListPTM.removeAll(toBeDeleteTaskList);
            taskListPTM.trimToSize();

            sqlConPTM.commit();

            int progress = countPtm * 20 / totalCountPtm;
            if (progress != lastProgressPtm) {
                logger.info("Post Searching {}%...", progress * 5);
                lastProgressPtm = progress;
            }

            if (countPtm == totalCountPtm) {
                break;
            }
            Thread.sleep(6000);
        }
        // shutdown threads.
        threadPoolPtm.shutdown();
        if (!threadPoolPtm.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolPtm.shutdownNow();
            if (!threadPoolPtm.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }

        sqlConPTM.commit();
        sqlConPTM.setAutoCommit(true);
        sqlConPTM.close();
        if (lockPtm.isLocked()) {
            lockPtm.unlock();
        }

        int fdrLevel = Integer.valueOf(parameterMap.get("fdr_level"));
        pfm(outputDir, allPeptideInfoMap, sqlPath, buildIndex.protSeqMap, massTool, hostName, varPtmCountMap, fdrLevel);

        logger.info("Saving results...");
    }

    private static void help() {
        String helpStr = "PIPI2: PIPI2: Sensitive Tag-based Database Search to Identify Peptides with Multiple PTMs.\r\n"
                + "Author: Shengzhi Lai\r\n"
                + "Email: slaiad@connect.ust.hk\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with PIPI2.\r\n"
                + "\t<data_file>: spectra data file\r\n"
                + "\texample: java -Xmx32g -jar PIPI2.jar parameter.def data.mgf \r\n";
        System.out.print(helpStr);
        System.exit(1);
    }

    class ScanRes {
        public int scanNum;
        public double qValue = -0.1;
        public List<CandiScore> peptideInfoScoreList;
        public double expMass;
        public int charge;
        public String scanName;
        ScanRes(String scanName, int scanNum, double expMass, List<CandiScore> peptideInfoScoreList, int charge){
            this.scanName = scanName;
            this.scanNum = scanNum;
            this.expMass = expMass;
            this.peptideInfoScoreList = peptideInfoScoreList;
            this.charge = charge;
        }
    }
    class CandiScore implements Comparable<CandiScore>{
        public String ptmContainingSeq;
        public PeptideInfo peptideInfo;
        public double pepScore;
        public double protScore = 0;
        public double varPtmTotalScore = 0;
        CandiScore(PeptideInfo peptideInfo, double pepScore, String ptmContainingSeq) {
            this.peptideInfo = peptideInfo;
            this.pepScore = pepScore;
            this.ptmContainingSeq = ptmContainingSeq;
        }

        public void setVarPtmTotalScore(Map<String, Double> varPtmStrScoreRefMap) { //varPtmScoreRefMap is like K(14.016)->1.2
            int startI = -1;
            double varPtmTotalScore = 0;
            int n_PtmOnPep = 0;
            for (int anyI = 0; anyI < ptmContainingSeq.length(); anyI++) {
                char thisChar = ptmContainingSeq.charAt(anyI);
                if (thisChar == '(') {
                    n_PtmOnPep++;
                    startI = anyI-1;
                } else if (thisChar == ')') {
                    String thisVarPtmStr = ptmContainingSeq.substring(startI, anyI + 1);
                    if (varPtmStrScoreRefMap.containsKey(thisVarPtmStr)) {
                        varPtmTotalScore += varPtmStrScoreRefMap.get(thisVarPtmStr);
                    }
                }
            }
            this.varPtmTotalScore = varPtmTotalScore/n_PtmOnPep;
        }
        public int compareTo(CandiScore o2) {
            if (this.protScore < o2.protScore) {
                return -1;
            } else if (this.protScore > o2.protScore) {
                return 1;
            } else {
                if (this.varPtmTotalScore < o2.varPtmTotalScore) {
                    return -1;
                } else if (this.varPtmTotalScore > o2.varPtmTotalScore) {
                    return 1;
                } else {
                    if (this.pepScore < o2.pepScore) {
                        return -1;
                    } else if (this.pepScore > o2.pepScore) {
                        return 1;
                    }
                }
            }
            return 0;
        }
    }
    private void pfm(String outputDir, Map<String, PeptideInfo> allPeptideInfoMap, String sqlPath, Map<String, String> protSeqMap, MassTool massTool, String hostName, Map<VarPtm, Integer> varPtmCountMap, int fdrLevel) throws IOException, SQLException {
        Map<String, Double> varPtmRefScoreMap = new HashMap<>();
        for (VarPtm varPtm : varPtmCountMap.keySet()) {
            if (varPtmCountMap.get(varPtm) == 1) {
                continue;
            }
            varPtmRefScoreMap.put(varPtm.site+"("+InferPTM.df3.format(varPtm.mass)+")", Math.sqrt(varPtmCountMap.get(varPtm)));
        }
        List<Map.Entry<String, Double>> testList = new ArrayList<>(varPtmRefScoreMap.entrySet());
        Collections.sort(testList, Map.Entry.comparingByValue(Comparator.reverseOrder()));
        int j = 0;
        for (Map.Entry<String, Double> entry : testList) {
            String varPtmStr = entry.getKey();
            if (j>maxNumVarPtmConsidered) {
                varPtmRefScoreMap.remove(varPtmStr);
            }
            j++;
        }

        //collect data
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, scanName, precursorCharge, precursorMass, peptide, theoMass, isDecoy, score, deltaLCn, deltaCn, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, peptideSet FROM spectraTable");
        Map<String, Map<String, Double>> protPepScoreMap = new HashMap<>();

        List<ScanRes> scanResList = new ArrayList<>();
        List<Double> topScoreList = new LinkedList<>();
        while (sqlResultSet.next()) {
            String topPeptide = sqlResultSet.getString("peptide");
            if (!sqlResultSet.wasNull()) {
                int charge = sqlResultSet.getInt("precursorCharge");
                double expMass = sqlResultSet.getDouble("precursorMass");
                double score = sqlResultSet.getDouble("score");
                String peptideSet = sqlResultSet.getString("peptideSet");
                int scanNum = sqlResultSet.getInt("scanNum");
                String scanName = sqlResultSet.getString("scanName");

                String[] candiSetStr = peptideSet.split(",");
                int numPep = candiSetStr.length/3;

                PeptideInfo pepInfo = allPeptideInfoMap.get(topPeptide.replaceAll("[^A-Z]+", ""));
                List<CandiScore> candiScoreList = new ArrayList<>();
                for (int i = 0; i < numPep; i++) {
                    String ptmContainingSeq = candiSetStr[3*i+0];
                    PeptideInfo candiPeptideInfo = allPeptideInfoMap.get(ptmContainingSeq.replaceAll("[^A-Z]+", ""));
                    double thisScore = Double.valueOf(candiSetStr[3*i+1]);
                    CandiScore candiScore = new CandiScore(candiPeptideInfo, thisScore, ptmContainingSeq);
                    candiScore.setVarPtmTotalScore(varPtmRefScoreMap);
                    candiScoreList.add(candiScore);
                }
                Collections.sort(candiScoreList, Comparator.comparing(o -> o.pepScore, Comparator.reverseOrder()));
                double topPepScore = candiScoreList.get(0).pepScore;
                topScoreList.add(topPepScore);
                Iterator<CandiScore> iter = candiScoreList.iterator();
                while (iter.hasNext()) {
                    CandiScore candiScore = iter.next();
                    if (candiScore.pepScore < 0.85 * topPepScore) {
                        iter.remove();
                    }
                }
                scanResList.add(new ScanRes(scanName, scanNum, expMass, candiScoreList, charge));

                for (String protId : pepInfo.protIdSet){
                    if (protPepScoreMap.containsKey(protId)){
                        Map<String, Double> pepScoreMap = protPepScoreMap.get(protId);
                        if (pepScoreMap.containsKey(topPeptide)){
                            pepScoreMap.put(topPeptide, Math.max(pepScoreMap.get(topPeptide), score));
                        } else {
                            pepScoreMap.put(topPeptide, score);
                        }
                    } else {
                        Map<String, Double> pepScoreMap = new HashMap<>();
                        pepScoreMap.put(topPeptide, score);
                        protPepScoreMap.put(protId, pepScoreMap);
                    }
                }
            }
        }
        sqlResultSet.close();
        sqlStatement.close();
        sqlConnection.close();

        Collections.sort(topScoreList);
        int numScores = topScoreList.size();
        double medianScore = numScores%2 == 0 ? (topScoreList.get(numScores/2) + topScoreList.get(numScores/2-1))/2 : topScoreList.get((numScores-1)/2);
        medianScore*=1.5;
        Map<String, Double> protScoreMap = new HashMap<>();
        for (String protId : protPepScoreMap.keySet()) {
            Map<String,Double> pepScoreMap = protPepScoreMap.get(protId);
            for (String pep : pepScoreMap.keySet()){
                if (pepScoreMap.get(pep) > medianScore) {
                    if (protScoreMap.containsKey(protId)){
                        protScoreMap.put(protId, protScoreMap.get(protId)+ pepScoreMap.get(pep));
                    } else {
                        protScoreMap.put(protId, pepScoreMap.get(pep));
                    }
                }
            }
        }
        for (String protId : protScoreMap.keySet()){
            protScoreMap.put(protId, protScoreMap.get(protId) / Math.log(protSeqMap.get(protId).length()));
        }
        for (ScanRes scanRes : scanResList) {
            List<CandiScore> candiScoreList = scanRes.peptideInfoScoreList;
            for (CandiScore candiScore : candiScoreList) {
                double protScoreForCand = -1;
                for (String protId : candiScore.peptideInfo.protIdSet) {
                    if (protScoreMap.getOrDefault(protId, 0.0) > protScoreForCand){
                        protScoreForCand = protScoreMap.getOrDefault(protId, 0.0);
                    }
                }
                candiScore.protScore = protScoreForCand;
            }
        }

        DecimalFormat df= new  DecimalFormat( ".00000" );

        //calculate ori FDR
        for (ScanRes scanRes : scanResList) {
            Collections.sort(scanRes.peptideInfoScoreList, Comparator.comparing(o -> o.pepScore, Comparator.reverseOrder()));
        }
        Collections.sort(scanResList, Comparator.comparing(o -> o.peptideInfoScoreList.get(0).pepScore, Comparator.reverseOrder()));
        double minQ = 1.0;
        boolean found = false;
        for (ScanRes scanRes : scanResList) {
            Collections.sort(scanRes.peptideInfoScoreList, Comparator.reverseOrder());
        }

        List<Double> fdrList = new ArrayList<>(scanResList.size());
        int numT = 0;
        int numD = 0;
        Collections.sort(scanResList, Comparator.comparing(o -> (o.peptideInfoScoreList.get(0).pepScore)*(o.peptideInfoScoreList.get(0).protScore+1), Comparator.reverseOrder())); // should still use peptideScore to do FDR
        for (ScanRes scanRes : scanResList) {
            CandiScore candiScore = scanRes.peptideInfoScoreList.get(0);
            if (candiScore.peptideInfo.isDecoy) {
                numD++;
            } else {
                numT++;
            }
            fdrList.add((double)numD/numT);
        }
        minQ = 1.0;
        found = false;
        for (int i = fdrList.size()-1; i >= 0; i--) {
            minQ = Math.min(fdrList.get(i), minQ);
            scanResList.get(i).qValue = minQ;
            if (!found && minQ < 0.01) {
                found = true;
            }
        }

        if (fdrLevel == 1) { // peptide level FDR
            List<ScanRes> copy_scanResList = new ArrayList<>(scanResList);
            Set<String> topPeps = new HashSet<>();
            List<Integer> indexesToDel = new ArrayList<>();

            for (int i = 0; i < copy_scanResList.size(); i++) {
                ScanRes thisRes = copy_scanResList.get(i);
                if (topPeps.contains(thisRes.peptideInfoScoreList.get(0).peptideInfo.freeSeq)) {
                    indexesToDel.add(i);
                } else {
                    topPeps.add(thisRes.peptideInfoScoreList.get(0).peptideInfo.freeSeq);
                }
            }
            indexesToDel.sort(Comparator.reverseOrder());
            for (int i : indexesToDel) {
                copy_scanResList.remove(i);
            }

            int newNumT = 0;
            int newNumD = 0;
            List<Double> newFdrList = new ArrayList<>();

            for (int i = 0; i < copy_scanResList.size(); i++) {
                ScanRes thisRes = copy_scanResList.get(i);
                if (thisRes.peptideInfoScoreList.get(0).peptideInfo.isDecoy) {
                    newNumD ++;
                } else {
                    newNumT++;
                }
                newFdrList.add(Math.min(1.0, (1+newNumD)/newNumT));
            }
            for (int i = newFdrList.size()-2; i >= 0; i--) {
                newFdrList.set(i, Math.min(newFdrList.get(i), newFdrList.get(i+1)));
            }

            Map<String, Double> pepQMap = new HashMap<>();
            for (int i = 0; i < copy_scanResList.size(); i++) {
                ScanRes thisRes = copy_scanResList.get(i);
                pepQMap.put(thisRes.peptideInfoScoreList.get(0).peptideInfo.freeSeq, newFdrList.get(i));
            }
            for (int i = 0; i < scanResList.size(); i++) {
                ScanRes thisRes = scanResList.get(i);
                thisRes.qValue = pepQMap.get(thisRes.peptideInfoScoreList.get(0).peptideInfo.freeSeq);
            }
        }

        List<Pair<Double, String>> finalExcelList = new ArrayList<>(scanResList.size());
        for (ScanRes scanRes : scanResList) {
            List<CandiScore> candiScoreList = scanRes.peptideInfoScoreList;
            CandiScore topCandi = candiScoreList.get(0);
            double theoMass = massTool.calResidueMass(topCandi.ptmContainingSeq) + massTool.H2O;
            double massDiff = getMassDiff(scanRes.expMass, theoMass, MassTool.C13_DIFF);
            double ppm = Math.abs(massDiff * 1e6 / theoMass);
            double finalScore = topCandi.pepScore*(topCandi.protScore+1);
            String finalStr = String.format(Locale.US, "%s,%d,%f,%d,%s,%s,%s,%s,%f,%f,%f,%d\n"
                    , scanRes.scanName, scanRes.scanNum, scanRes.qValue, topCandi.peptideInfo.isDecoy ? 0 : 1, df.format(finalScore), topCandi.ptmContainingSeq, topCandi.peptideInfo.freeSeq
                    , String.join(";",topCandi.peptideInfo.protIdSet), ppm, theoMass, scanRes.expMass, scanRes.charge
            );

            finalExcelList.add(new Pair(finalScore, finalStr));
        }
        Collections.sort(finalExcelList, Comparator.comparing(o -> o.getFirst(), Comparator.reverseOrder()));
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputDir+"PIPI."+hostName+".csv"));

        writer.write("scanName,scanNum,qValue,TorD,score,modifiedSeq,freeSeq,proteins,ppm,theoMass,expMass,charge\n");
        for (Pair<Double, String> pair : finalExcelList) {
            writer.write(pair.getSecond());
        }
        writer.close();
    }

    public static double getMassDiff(double expMass, double theoMass, double C13Diff) {
        double massDiff1 = expMass - theoMass;
        double massDiff2 = expMass - theoMass - C13Diff;
        double massDiff3 = expMass - theoMass - 2 * C13Diff;
        double absMassDiff1 = Math.abs(massDiff1);
        double absMassDiff2 = Math.abs(massDiff2);
        double absMassDiff3 = Math.abs(massDiff3);

        if ((absMassDiff1 <= absMassDiff2) && (absMassDiff1 <= absMassDiff2)) {
            return massDiff1;
        } else if ((absMassDiff2 <= absMassDiff1) && (absMassDiff2 <= absMassDiff3)) {
            return massDiff2;
        } else {
            return massDiff3;
        }
    }
}
