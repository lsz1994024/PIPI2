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

import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.SpecProcessor;
import ProteomicsLibrary.Types.SparseBooleanVector;
import ProteomicsLibrary.Types.SparseVector;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.FM.FMIndex;
import proteomics.FM.FMRes;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

import static proteomics.PIPI.*;
import static proteomics.PTM.InferPTM.*;
import static proteomics.PTM.InferPTM.N_PART;
import static proteomics.Segment.InferSegment.*;

public class PreSearch implements Callable<PreSearch.Entry> {
    private static final Logger logger = LoggerFactory.getLogger(PreSearch.class);

    private static final int candisNum = 20;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms1Tolerance;
    private final double leftInverseMs1Tolerance;
    private final double rightInverseMs1Tolerance;
    private final int ms1ToleranceUnit;
    private final double minPtmMass;
    private final double maxPtmMass;
    private final int localMaxMs2Charge;
    private final JMzReader spectraParser;
    private final double minClear;
    private final double maxClear;
    private final ReentrantLock lock;
    private final String scanName;
    private final int precursorCharge;
    private final double precursorMass;
    private final SpecProcessor specProcessor;
    private final int scanNum;
    private String truth;
    private final double ms2Tolerance;
    private final InferPTM inferPTM;
    private final int minPepLen;
    private final int maxPepLen;
    public PreSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms2Tolerance, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance
            , int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int localMaxMs2Charge
            , JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanName, int precursorCharge, double precursorMass
            , SpecProcessor specProcessor, String truth, int minPepLen, int maxPepLen) {

        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms1Tolerance = ms1Tolerance;
        this.leftInverseMs1Tolerance = leftInverseMs1Tolerance;
        this.rightInverseMs1Tolerance = rightInverseMs1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.localMaxMs2Charge = localMaxMs2Charge;
        this.spectraParser = spectraParser;
        this.minClear = minClear;
        this.maxClear = maxClear;
        this.lock = lock;
        this.scanName = scanName;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.specProcessor = specProcessor;
        this.scanNum = scanNum;
        this.truth = truth;
        this.ms2Tolerance = ms2Tolerance;
        this.inferPTM = buildIndex.getInferPTM();
        this.minPepLen = minPepLen;
        this.maxPepLen = maxPepLen;
    }

    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            rawPLMap = spectraParser.getSpectrumById(scanName.split("\\.")[1]).getPeakList();
        } finally {
            lock.unlock();
        }
        double ms1TolAbs = Double.parseDouble(InferPTM.df3.format(precursorMass*ms1Tolerance/1000000));
        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN, ms2Tolerance);

        if (plMap.isEmpty()) return null;
        // Coding
        SparseVector expProcessedPL;
        if (PIPI.useXcorr) {
            expProcessedPL = specProcessor.prepareXCorr(plMap, false);
        } else {
            expProcessedPL = specProcessor.digitizePL(plMap);
        }
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);

        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, minTagLenToExtract,maxTagLenToExtract);

        List<ExpTag> cleanedAllLongTagList = inferSegment.cleanAbundantTagsPrefix(allLongTagList, minTagLenToExtract);

        allLongTagList = cleanedAllLongTagList;

        if (allLongTagList.isEmpty())  return null;

        double totalMass = precursorMass + 2 * MassTool.PROTON;

        FMIndex fmIndex = buildIndex.fmIndexReduced;

        Map<String, List<Peptide>> resPeptideListMap = new HashMap<>();
        Map<String, PeptideInfo> peptideInfoMap = new HashMap<>(50000);

        Set<String> searchedTagStrSet = new HashSet<>();
        int minTagLen = 4;

        int n_tags = 0;
        for (ExpTag tagInfo : allLongTagList.subList(0, Math.min(10, allLongTagList.size()))){

            minTagLen = tagInfo.size() > 4 ? 5 : 4;
            if (buildIndex.posProtMapReduced.size() < 5000) {
                minTagLen = 3;
            }
            String tagStr = tagInfo.getFreeAaString();
            String revTagStr = new StringBuilder(tagStr).reverse().toString();

            if (tagInfo.isNorC == N_TAG) { //n tag
                String tagStrMzStr = tagStr + df3.format(tagInfo.getHeadLocation());
                if (!searchedTagStrSet.contains(tagStrMzStr)) {
                    char[] tagChar = tagStr.toCharArray();

                    Set<String> protIdSetByThisTag = new HashSet<>();
                    int n_res = searchAndSaveFuzzy(scanNum, tagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, tagChar, minTagLen, expProcessedPL,finalPlMap, true);
                    searchedTagStrSet.add(tagStrMzStr);
                    if (n_res > 100) {
                        continue;
                    }
                    if (tagStr.length() > minTagLen+1){ // if the tag was already long i.e. there is space to sub
                        for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min(2,tagChar.length-minTagLen); n_cTermAaToCut++){
                            String subTagStr = tagStr.substring(0, tagStr.length()-n_cTermAaToCut);
                            char[] subTagChar = subTagStr.toCharArray();
                            ExpTag subTagInfo = tagInfo.subTag(0,tagStr.length()-n_cTermAaToCut);
                            if (!searchedTagStrSet.contains(tagStr)) {
                                int numResSub = searchAndSaveFuzzy(scanNum ,subTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, subTagChar, minTagLen, expProcessedPL, finalPlMap, false);
                                searchedTagStrSet.add(tagStr);
                                if (numResSub > 100) {
                                    break;
                                }
                            }
                        }
                    }
                }
            } else if (tagInfo.isNorC == C_TAG) { // c tag
                char[] revTagChar = revTagStr.toCharArray();
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                String revTagStrMzStr = revTagStr + df3.format(revTagInfo.getHeadLocation());
                if (!searchedTagStrSet.contains(revTagStrMzStr)) {
                    Set<String> protIdSetByThisTag = new HashSet<>();
                    int n_res = searchAndSaveFuzzy(scanNum ,revTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, revTagChar, minTagLen, expProcessedPL, finalPlMap, true);
                    searchedTagStrSet.add(revTagStrMzStr);
                    if (n_res > 100) {
                        continue;
                    }
                    if (revTagStr.length() > minTagLen+1){ // if the tag was already long i.e. there is space to sub
                        for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min(2,revTagChar.length-minTagLen); n_cTermAaToCut++){
                            String subRevTagStr = revTagStr.substring(0, tagStr.length()-n_cTermAaToCut);
                            char[] subRevTagChar = subRevTagStr.toCharArray();
                            ExpTag subRevTagInfo = revTagInfo.subTag(0,tagStr.length()-n_cTermAaToCut);
                            if (!searchedTagStrSet.contains(subRevTagStr)) {
                                subRevTagInfo.isNorC = NON_NC_TAG;
                                int numResSub = searchAndSaveFuzzy(scanNum ,subRevTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, subRevTagChar, minTagLen, expProcessedPL, finalPlMap, false);
                                searchedTagStrSet.add(subRevTagStr);
                                if (numResSub > 100) {
                                    break;
                                }
                            }
                        }
                    }
                }
            } else { // non-nc tag
                char[] tagChar = tagStr.toCharArray();
                String tagStrMzStr = tagStr + df3.format(tagInfo.getHeadLocation());
                if (!searchedTagStrSet.contains(tagStrMzStr)) {
                    Set<String> protIdSetByThisTag = new HashSet<>();
                    int n_res = searchAndSaveFuzzy(scanNum ,tagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, tagChar, minTagLen, expProcessedPL,finalPlMap, true);
                    searchedTagStrSet.add(tagStrMzStr);
                    if (n_res < 100 ) {
                        if (tagStr.length() > minTagLen+1){
                            for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min(2,tagChar.length-minTagLen); n_cTermAaToCut++){
                                String subTagStr = tagStr.substring(0, tagStr.length()-n_cTermAaToCut);
                                ExpTag subTagInfo = tagInfo.subTag(0,tagStr.length()-n_cTermAaToCut);
                                char[] subTagChar = subTagStr.toCharArray();
                                int numResSub1 = 0;
                                if (!searchedTagStrSet.contains(subTagStr)) {
                                    numResSub1 = searchAndSaveFuzzy(scanNum ,subTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, subTagChar, minTagLen, expProcessedPL,finalPlMap, false);
                                    searchedTagStrSet.add(subTagStr);
                                    if (numResSub1 > 100) {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                char[] revTagChar = revTagStr.toCharArray();
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                String revTagStrMzStr = revTagStr + df3.format(revTagInfo.getHeadLocation());
                if (!searchedTagStrSet.contains(revTagStrMzStr)) {
                    Set<String> protIdSetByThisTag = new HashSet<>();
                    int n_res = searchAndSaveFuzzy(scanNum ,revTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, revTagChar, minTagLen, expProcessedPL, finalPlMap, true);
                    searchedTagStrSet.add(revTagStrMzStr);
                    if (n_res < 100 ) {
                        if (revTagStr.length() > minTagLen+1){
                            for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min(2,revTagChar.length-minTagLen); n_cTermAaToCut++){
                                String subTagStr = revTagStr.substring(0, revTagStr.length()-n_cTermAaToCut);
                                ExpTag subTagInfo = tagInfo.revTag(totalMass).subTag(0,revTagStr.length()-n_cTermAaToCut);

                                //sub forward
                                char[] subTagChar = subTagStr.toCharArray();
                                int numResSub1 = 0;
                                if (!searchedTagStrSet.contains(subTagStr)) {
                                    numResSub1 = searchAndSaveFuzzy(scanNum ,subTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, subTagChar, minTagLen, expProcessedPL, finalPlMap, false);
                                    searchedTagStrSet.add(subTagStr);
                                    if (numResSub1 > 100) {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            n_tags++;
        }

        TreeSet<Peptide> peptideSet = new TreeSet<>(Collections.reverseOrder());
        if (!resPeptideListMap.isEmpty()) {
            for (String pepSeq : resPeptideListMap.keySet()) {
                List<Peptide> peptideList = resPeptideListMap.get(pepSeq);
                Collections.sort(peptideList, Comparator.comparing(o->o.getScore(), Comparator.reverseOrder()));
                double topScore = peptideList.get(0).getScore();
                for (Peptide candiPep : peptideList) {
                    if (candiPep.getScore() > 0* topScore){
                        if (candiPep.getScore() > 0) {
                            if (peptideSet.size() < candisNum) {
                                peptideSet.add(candiPep);
                            } else if (candiPep.compareTo( peptideSet.last() ) == 1) {
                                peptideSet.pollLast();
                                peptideSet.add(candiPep);
                            }
                        }
                    }
                }
            }
        }
        if (peptideSet.isEmpty()) {
            return null;
        }
        List<Peptide> pepList = new ArrayList<>(peptideSet);
        Set<Integer> pepIdsToRemove = new HashSet<>();
        for (int j = pepList.size()-1; j >= 1; j--) {
            if (pepIdsToRemove.contains(j)) continue;
            for (int i = 0; i < j; i++) {
                if (pepIdsToRemove.contains(i)) continue;
                if (isHomo(pepList.get(i), pepList.get(j), peptideInfoMap) || pepList.get(j).getScore() > 0.6*pepList.get(i).getScore()) {
                    int iPriority = pepList.get(i).getPriority();
                    int jPriority = pepList.get(j).getPriority();
                    if (iPriority < jPriority) {
                        pepIdsToRemove.add(i);
                    } else if (iPriority > jPriority) {
                        pepIdsToRemove.add(j);
                    } else {// iPriority == jPriority

                        if (onlyDifferUnsettledPtm(pepList.get(i), pepList.get(j))) {
                            pepIdsToRemove.add(j);
                        }
                    }
                }
            }
            if (pepList.get(j).getPriority() < 0 && isPtmSimuTest ) {
                pepIdsToRemove.add(j);
            }
        }
        if (pepList.get(0).getPriority() < 0 && isPtmSimuTest ) {
            pepIdsToRemove.add(0);
        }
        List<Peptide> newPepList = new ArrayList<>();
        for (int id = 0; id < pepList.size(); id++){
            if (pepIdsToRemove.contains(id)) continue;
            newPepList.add(pepList.get(id));
        }
        if (newPepList.isEmpty()) {
            return null;
        }

        Peptide[] peptideArray = newPepList.toArray(new Peptide[0]);
        Peptide topPep = peptideArray[0];

        Entry entry;
        String pepSetString = "";
        for (Peptide peptide : peptideArray){
            PeptideInfo peptideInfo = peptideInfoMap.get(peptide.getFreeSeq());
            pepSetString += peptide.getVarPtmContainingSeqNow() + "," + peptide.getScore() + "," + String.join("_", peptideInfo.protIdSet) +",";
        }

        boolean shouldPtm = Math.abs(precursorMass-massTool.calResidueMass(topPep.getFreeSeq()) - massTool.H2O) > 0.01;
        boolean hasPTM = topPep.hasVarPTM();
        int ptmNum = 0;
        boolean isSettled = true;
        double totalPtmMass = 0;
        if (hasPTM) {
            ptmNum = topPep.getVarPTMs().size();
            for (double mass : topPep.getVarPTMs().values()){
                totalPtmMass += mass;
            }
        }
        isSettled = Math.abs(totalPtmMass-(precursorMass-massTool.calResidueMass(topPep.getFreeSeq()) - massTool.H2O)) <= 0.01;

        double deltaLCn = 1;
        if (peptideArray.length > candisNum - 1) {
            deltaLCn = (peptideArray[0].getScore() - peptideArray[candisNum - 1].getScore()) / peptideArray[0].getScore();
        }
        double deltaCn = 1;
        if (peptideArray.length > 1) {
            for(int i = 0; i < peptideArray.length; i++) {
                if (peptideArray[i].getScore() != peptideArray[0].getScore()){
                    deltaCn = (peptideArray[0].getScore() - peptideArray[i].getScore()) / peptideArray[0].getScore();
                    break;
                }
            }
        }
        String otherPtmPatterns = "-";
        entry = new PreSearch.Entry(
                scanNum, scanName, shouldPtm ? 1 : 0, hasPTM ? 1 : 0, ptmNum, isSettled ? 1 : 0
                , precursorCharge, precursorMass, buildIndex.getLabelling(), topPep.getPtmContainingSeq(buildIndex.returnFixModMap())
                , topPep.getTheoMass(), topPep.isDecoy() ? 1 : 0, topPep.getScore(), deltaLCn, deltaCn
                , topPep.getMatchedPeakNum(), topPep.getIonFrac(), topPep.getMatchedHighestIntensityFrac()
                , topPep.getExplainedAaFrac(), otherPtmPatterns, topPep.getaScore(), ""
                , pepSetString.substring(0, pepSetString.length()-1)
        );
        entry.topPeptide = topPep;
        entry.topPeptide.precursorMass = precursorMass;
        entry.varPtmList.addAll(topPep.posVarPtmResMap.values());
        for (Peptide peptide : peptideArray) {
            entry.peptideInfoMapForRef.put(peptide.getFreeSeq(), peptideInfoMap.get(peptide.getFreeSeq()));
        }

        return entry;
    }
    private boolean isHomo(Peptide p1, Peptide p2, Map<String, PeptideInfo> peptideInfoMap) {

        Set<String> temp1 = peptideInfoMap.get(p1.getFreeSeq()).protIdSet;
        Set<String> temp2 = peptideInfoMap.get(p2.getFreeSeq()).protIdSet;
        Set<String> set = new HashSet<>(temp1);
        set.retainAll(temp2);
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getFreeSeq());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getFreeSeq());
        return sbv1.dot(sbv2) > 0.3*Math.min(p1.getFreeSeq().length(), p2.getFreeSeq().length());
    }

    private boolean onlyDifferUnsettledPtm(Peptide p1, Peptide p2) {
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getFreeSeq());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getFreeSeq());
        if (sbv1.dot(sbv2) < 0.3*Math.min(p1.getFreeSeq().length(), p2.getFreeSeq().length())){
            return false;
        }

        if (p1.posVarPtmResMap.size() != p2.posVarPtmResMap.size()) {
            return false;
        }
        if (!p1.posVarPtmResMap.keySet().containsAll(p2.posVarPtmResMap.keySet())){
            return false;
        }
        byte n_SameMass = 0;
        for (int pos : p2.posVarPtmResMap.keySet()) {
            if (p1.posVarPtmResMap.get(pos).mass == p2.posVarPtmResMap.get(pos).mass) {
                n_SameMass++;
            } else if (Math.abs(p1.posVarPtmResMap.get(pos).mass - p2.posVarPtmResMap.get(pos).mass) < 0.02) {
                if (p1.posVarPtmResMap.get(pos).classification.contains("PIPI_unsettled") || p2.posVarPtmResMap.get(pos).classification.contains("PIPI_unsettled")) {
                    n_SameMass++;
                }
            }
        }
        return n_SameMass == p1.posVarPtmResMap.size();
    }


    private int searchAndSaveFuzzy(int scanNum, ExpTag tagInfo, double ms1TolAbs, Map<String, List<Peptide>> resPeptideListMap, Map<String, PeptideInfo> peptideInfoMap
            , Set<String> protIdSetByThisTag, FMIndex fmIndex, char[] tagChar, int minTagLen, SparseVector expProcessedPL,TreeMap<Double, Double> plMap, boolean isFirstTime) throws CloneNotSupportedException {

        int numRes = 0;
        FMRes searchRes = fmIndex.fmSearchFuzzy(tagChar);
        if (searchRes == null) {
            return 0;
        }
        numRes = searchRes.ep-searchRes.sp+1;
        int solCount = 0;
        if (searchRes.settled) {
            for (int ii = searchRes.sp; ii <= searchRes.ep; ii++) {
                if (solCount > 100) break;
                int absTagPos = fmIndex.SA[ii];
                int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);

                String protId = buildIndex.posProtMapReduced.get(dotIndex);
                if (protIdSetByThisTag.contains(protId)){
                    continue;
                }
                solCount++;
                protIdSetByThisTag.add(protId);
                int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;
                updateCandiList(scanNum, protId, relPos, tagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, expProcessedPL, plMap);
            }
        } else {
            if (!isFirstTime) return 0;
            int matchedPos = searchRes.matchedPos;
            if (tagInfo.size()-matchedPos < minTagLen) return 0;
            for (int ii = searchRes.sp; ii <= searchRes.ep; ii++) {
                if (solCount > 100) break;
                int absTagPos = fmIndex.SA[ii];
                int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);
                String protId = buildIndex.posProtMapReduced.get(dotIndex);
                String protSeq = buildIndex.protSeqMap.get(protId);
                int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;
                solCount++;
                updateCandiList(scanNum, protId, relPos, tagInfo.subTag(matchedPos, tagChar.length), ms1TolAbs, resPeptideListMap, peptideInfoMap, expProcessedPL, plMap);
            }

        }
        return solCount;
    }

    private void updateCandiList(int scanNum, String protId, int tagPosInProt, ExpTag finderTag, double ms1TolAbs, Map<String, List<Peptide>> resPeptideListMap
            , Map<String, PeptideInfo> peptideInfoMap, SparseVector expProcessedPL, TreeMap<Double, Double> plMap) throws CloneNotSupportedException {

        double tagCMass = finderTag.getTailLocation() + massTool.H2O-MassTool.PROTON;
        String protSeq = buildIndex.protSeqMap.get(protId);
        Map<Integer, Double> cPoscMassMap = new HashMap<>();
        if ( (tagPosInProt+finderTag.size() == protSeq.length() && finderTag.isNorC != C_TAG)
                || (tagPosInProt == 0 && finderTag.isNorC != N_TAG)
                || tagPosInProt < 0
                || tagPosInProt+finderTag.size() > protSeq.length()
        ) {
            return;
        }

        int missCleav = getNumOfMissCleavSite(finderTag.getFreeAaString());
        if (finderTag.isNorC == C_TAG || Math.abs(tagCMass - precursorMass) < ms1TolAbs){
            if (isKR(protSeq.charAt(tagPosInProt+finderTag.size()-1)) && Math.abs(tagCMass - precursorMass) < ms1TolAbs) {
                cPoscMassMap.put(tagPosInProt+finderTag.size()-1, tagCMass);
            }
        } else {
            for (int i = tagPosInProt+finderTag.size(); i < protSeq.length(); i++) {  //
                if (isX(protSeq.charAt(i))) break;
                tagCMass += massTool.getMassTable().get(protSeq.charAt(i));
                if (tagCMass > precursorMass+maxPtmMass) break;

                if (isKR(protSeq.charAt(i))) {
                    missCleav++;
                }
                if (tagCMass >= precursorMass-maxPtmMass && isKR(protSeq.charAt(i))) {
                    cPoscMassMap.put(i, tagCMass);
                }
                if (missCleav > massTool.missedCleavage  && !isPtmSimuTest) {
                    break;
                }
            }
        }

        for (int cPos : cPoscMassMap.keySet()) {
            double cDeltaMass = precursorMass - cPoscMassMap.get(cPos);
            List<TreeMap<Integer, VarPtm>> cPosVarPtmResMap_List = new LinkedList<>();
            char rightFlank;
            byte isProtNorC_Term = NON_TERM_PROT;
            if (cPos == protSeq.length()-1) {
                rightFlank = '-';
                isProtNorC_Term = C_TERM_PROT;
            } else {
                rightFlank = protSeq.charAt(cPos+1);
            }
            String cPartSeq = protSeq.substring(tagPosInProt+finderTag.size(), cPos+1);
            double cCutMass = finderTag.getTailLocation() - MassTool.PROTON;
            boolean isCTermFree = false;
            if (Math.abs(cDeltaMass) > 0.1) {
                Set<Integer> fixModIdxes = inferPTM.getFixModIdxes(cPartSeq);
                Map<Integer, VarPtm[]> partPosVarModArrayMap = inferPTM.getIdxVarModMapNew(cPartSeq, fixModIdxes, C_PART, isProtNorC_Term);
                ModPepPool cPartPepsWithPtm = inferPTM.settlePtmOnSide(scanNum, expProcessedPL, plMap, precursorMass, cPartSeq, false,
                        partPosVarModArrayMap, cCutMass, cDeltaMass, precursorCharge, C_PART,ms1TolAbs);

                if (cPartPepsWithPtm.peptideTreeSet.isEmpty()) {
                    continue;
                }
                double topScore = cPartPepsWithPtm.getTopPepPtn().getScore();
                boolean shouldSkip = true;
                for (Peptide peptide : cPartPepsWithPtm.peptideTreeSet){
                    if (peptide.getScore() >= 0.5 *topScore){
                        if (!peptide.posVarPtmResMap.isEmpty()){
                            cPosVarPtmResMap_List.add(peptide.posVarPtmResMap);
                            shouldSkip = false;
                        }
                    }
                }
                if (shouldSkip) {
                    continue;
                }
            } else {
                isCTermFree = true;
            }


            double nDeltaMass = finderTag.getHeadLocation() -MassTool.PROTON;
            missCleav = getNumOfMissCleavSite(protSeq.substring(tagPosInProt, cPos)) ;

            int min_nPos = (finderTag.isNorC == N_TAG && Math.abs(finderTag.getHeadLocation()-MassTool.PROTON) <= 0.02 ) ? tagPosInProt : 0;
            int max_nPos = (finderTag.isNorC == N_TAG && Math.abs(finderTag.getHeadLocation()-MassTool.PROTON) <= 0.02 ) ? tagPosInProt : tagPosInProt-1;
            for (int nPos = max_nPos; nPos >= min_nPos; nPos--) {
                if (cPos+1-nPos < minPepLen){
                    continue;
                }
                else if ( cPos+1-nPos > maxPepLen) {
                    break;
                }
                if (isX(protSeq.charAt(nPos))) break;
                if (nPos < tagPosInProt) {
                    nDeltaMass -= massTool.getMassTable().get(protSeq.charAt(nPos));
                }
                if (nDeltaMass < minPtmMass) break;
                if (nDeltaMass > maxPtmMass)
                    continue;

                if (nTermSpecific) {
                    if (nPos != 0 && !isKR(protSeq.charAt(nPos - 1))) {
                        continue;
                    }
                }

                missCleav = getNumOfMissCleavSite(protSeq.substring(nPos, cPos));
                if (missCleav > massTool.missedCleavage && !isPtmSimuTest) {
                    break;
                }
                if (cPos + 1 - nPos < 6) {
                    continue;
                }

                String nPartSeq = protSeq.substring(nPos, tagPosInProt);
                char leftFlank;
                isProtNorC_Term = NON_TERM_PROT;
                if (nPos == 0 || (nPos == 1 && protSeq.charAt(0) == 'M')) {
                    leftFlank = '-';
                    isProtNorC_Term = N_TERM_PROT;
                } else {
                    leftFlank = protSeq.charAt(nPos - 1);
                }
                double nCutMass = precursorMass+MassTool.PROTON - finderTag.getHeadLocation() - massTool.H2O;

                List<TreeMap<Integer, VarPtm>> nPosVarPtmResMap_List = new LinkedList<>();

                boolean isNTermFree = false;
                if (Math.abs(nDeltaMass) > 0.1) {
                    Set<Integer> fixModIdxes = inferPTM.getFixModIdxes(nPartSeq);
                    Map<Integer, VarPtm[]> partPosVarModArrayMap = inferPTM.getIdxVarModMapNew(nPartSeq, fixModIdxes, N_PART, isProtNorC_Term);
                    ModPepPool nPartpepsWithPtm = inferPTM.settlePtmOnSide(scanNum, expProcessedPL, plMap, precursorMass, nPartSeq, false,
                            partPosVarModArrayMap, nCutMass, nDeltaMass, precursorCharge, N_PART, ms1TolAbs);

                    if (nPartpepsWithPtm.peptideTreeSet.isEmpty()) {
                        continue;
                    }
                    double topScore = nPartpepsWithPtm.getTopPepPtn().getScore();
                    boolean shouldSkip = true;
                    for (Peptide peptide : nPartpepsWithPtm.peptideTreeSet) {
                        if (peptide.getScore() >= 0.5 * topScore) {
                            if (!peptide.posVarPtmResMap.isEmpty()) {
                                nPosVarPtmResMap_List.add(peptide.posVarPtmResMap);
                                shouldSkip = false;
                            }
                        }
                    }
                    if (shouldSkip) {
                        continue;
                    }
                } else {
                    isNTermFree = true;
                }
                int tagPosInPep = tagPosInProt - nPos;
                String freePepSeq = protSeq.substring(nPos, cPos + 1);
                StringBuilder ptmPepSeqSB = new StringBuilder(freePepSeq);
                ptmPepSeqSB.replace(tagPosInPep, tagPosInPep + finderTag.size(), finderTag.getPtmAaString());
                List<PosMassMap> posMassMap_List = new LinkedList<>();
                List<TreeMap<Integer, VarPtm>> posVarPtmResMap_List = new LinkedList<>();
                if (!cPosVarPtmResMap_List.isEmpty() && !nPosVarPtmResMap_List.isEmpty()) {
                    for (TreeMap<Integer, VarPtm> cPosVarPtmResMap : cPosVarPtmResMap_List){
                        for (TreeMap<Integer, VarPtm> nPosVarPtmResMap : nPosVarPtmResMap_List) {
                            PosMassMap fullPosMassMap = new PosMassMap(freePepSeq.length());
                            TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();
                            for (int pos : nPosVarPtmResMap.keySet()) {
                                fullPosMassMap.put(pos, nPosVarPtmResMap.get(pos).mass);
                                posVarPtmResMap.put(pos, nPosVarPtmResMap.get(pos));
                            }
                            for (int pos : cPosVarPtmResMap.keySet()) {
                                fullPosMassMap.put(pos + tagPosInPep + finderTag.size(), cPosVarPtmResMap.get(pos).mass);
                                posVarPtmResMap.put(pos + tagPosInPep + finderTag.size(), cPosVarPtmResMap.get(pos));
                            }
                            int idOfAa = -1;
                            for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
                                if (Character.isUpperCase(aaChar)) {
                                    idOfAa += 1;
                                } else {
                                    fullPosMassMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
                                    posVarPtmResMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar));
                                }
                            }
                            posVarPtmResMap_List.add(posVarPtmResMap);
                            posMassMap_List.add(fullPosMassMap);
                        }
                    }
                } else if (cPosVarPtmResMap_List.isEmpty() && !nPosVarPtmResMap_List.isEmpty()) {
                    for (TreeMap<Integer, VarPtm> nPosVarPtmResMap : nPosVarPtmResMap_List) {
                        PosMassMap fullPosMassMap = new PosMassMap(freePepSeq.length());
                        TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();
                        for (int pos : nPosVarPtmResMap.keySet()) {
                            fullPosMassMap.put(pos, nPosVarPtmResMap.get(pos).mass);
                            posVarPtmResMap.put(pos, nPosVarPtmResMap.get(pos));
                        }
                        int idOfAa = -1;
                        for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
                            if (Character.isUpperCase(aaChar)) {
                                idOfAa += 1;
                            } else {
                                fullPosMassMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
                                posVarPtmResMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar));
                            }
                        }
                        posVarPtmResMap_List.add(posVarPtmResMap);
                        posMassMap_List.add(fullPosMassMap);
                    }
                }else if (!cPosVarPtmResMap_List.isEmpty() && nPosVarPtmResMap_List.isEmpty()){
                    for (TreeMap<Integer, VarPtm> cPosVarPtmResMap : cPosVarPtmResMap_List){
                        PosMassMap fullPosMassMap = new PosMassMap(freePepSeq.length());
                        TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();
                        for (int pos : cPosVarPtmResMap.keySet()) {
                            fullPosMassMap.put(pos + tagPosInPep + finderTag.size(), cPosVarPtmResMap.get(pos).mass);
                            posVarPtmResMap.put(pos + tagPosInPep + finderTag.size(), cPosVarPtmResMap.get(pos));
                        }
                        int idOfAa = -1;
                        for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
                            if (Character.isUpperCase(aaChar)) {
                                idOfAa += 1;
                            } else {
                                fullPosMassMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
                                posVarPtmResMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar));
                            }
                        }
                        posVarPtmResMap_List.add(posVarPtmResMap);
                        posMassMap_List.add(fullPosMassMap);
                    }
                }

                if (posVarPtmResMap_List.isEmpty()){
                    Peptide peptide = new Peptide(freePepSeq, true, massTool);
                    peptide.tagPosInPep = tagPosInPep;
                    peptide.ptmSeq = ptmPepSeqSB.toString();
                    peptide.finderTag = finderTag;
                    peptide.cDeltaMass = cDeltaMass;
                    peptide.nDeltaMass = nDeltaMass;

                    peptide.setScore(massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL, peptide.matchedBions, peptide.matchedYions));//todo decide the penalty

                    List<Peptide> peptideList = resPeptideListMap.get(freePepSeq);
                    if (peptideList == null) {
                        peptideList = new ArrayList<>();
                        peptideList.add(peptide);
                        resPeptideListMap.put(freePepSeq, peptideList);
                    } else {
                        peptideList.add(peptide);
                    }
                }else {
                    for (int i = 0; i < posVarPtmResMap_List.size(); i++) {
                        Peptide peptide = new Peptide(freePepSeq, true, massTool);
                        peptide.tagPosInPep = tagPosInPep;
                        peptide.ptmSeq = ptmPepSeqSB.toString();
                        peptide.finderTag = finderTag;
                        peptide.cDeltaMass = cDeltaMass;
                        peptide.nDeltaMass = nDeltaMass;

                        PosMassMap fullPosMassMap = posMassMap_List.get(i);
                        TreeMap<Integer, VarPtm> posVarPtmResMap = posVarPtmResMap_List.get(i);
                        if (!fullPosMassMap.isEmpty()) {
                            peptide.setVarPTM(fullPosMassMap);
                            peptide.posVarPtmResMap = posVarPtmResMap;
                        }
                        peptide.setScore(massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL, peptide.matchedBions, peptide.matchedYions));//todo decide the penalty

                        List<Peptide> peptideList = resPeptideListMap.get(freePepSeq);
                        if (peptideList == null) {
                            peptideList = new ArrayList<>();
                            peptideList.add(peptide);
                            resPeptideListMap.put(freePepSeq, peptideList);
                        } else {
                            peptideList.add(peptide);
                        }
                    }
                }
                if (peptideInfoMap.containsKey(freePepSeq)) {
                    PeptideInfo peptideInfo = peptideInfoMap.get(freePepSeq);
                    if (!peptideInfo.protIdSet.contains(protId)) {
                        if (peptideInfo.leftFlank != '-' && peptideInfo.rightFlank != '-') {
                            if (rightFlank == '-' || leftFlank == '-') {
                                peptideInfo.leftFlank = leftFlank;
                                peptideInfo.rightFlank = rightFlank;
                            }
                        }
                        peptideInfo.protIdSet.add(protId);
                        if (!protId.startsWith("DECOY_")) {
                            peptideInfo.isDecoy = false;
                        }
                    }
                } else {
                    PeptideInfo peptideInfo = new PeptideInfo(freePepSeq, protId.startsWith("DECOY_"), leftFlank, rightFlank);
                    peptideInfo.protIdSet.add(protId);
                    peptideInfoMap.put(freePepSeq, peptideInfo);
                }

            }
        }

        return;
    }

    private int getNumOfMissCleavSite(String seq) {
        String str1 = seq.substring(0,seq.length());
        String str2 = str1.replaceAll("[KR]","");
        return str1.length()-str2.length();
    }
    private boolean isKR(char aa){
        return aa == 'K' || aa == 'R';
    }

    private boolean isX(char aa){
        return aa == 'X';
    }
    public class Entry {
        public Map<String, PeptideInfo> peptideInfoMapForRef = new HashMap<>();
        public Peptide topPeptide = null;
        final int scanNum;
        final String scanName;
        final int precursorCharge;
        final double precursorMass;
        final String labelling;
        public final String peptide;
        final double theoMass;
        final int isDecoy;
        public final double score;
        final double deltaLCn;
        final double deltaCn;
        final int matchedPeakNum;
        final double ionFrac;
        final double matchedHighestIntensityFrac;
        final double explainedAaFrac;
        final String otherPtmPatterns;
        final String aScore;
        final String candidates;
        final String peptideSet;
        final int hasPTM;
        final int ptmNum;
        final int isSettled;
        final int shouldPtm;

        List<VarPtm> varPtmList = new ArrayList<>();
        Entry(int scanNum, String scanName, int shouldPtm, int hasPTM, int ptmNum, int isSetteld, int precursorCharge, double precursorMass
                ,String labelling, String peptide, double theoMass, int isDecoy
                , double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac
                , double explainedAaFrac, String otherPtmPatterns, String aScore, String candidates, String peptideSet) {
            this.scanNum = scanNum;
            this.scanName = scanName;
            this.shouldPtm = shouldPtm;
            this.hasPTM = hasPTM;
            this.ptmNum = ptmNum;
            this.isSettled = isSetteld;
            this.precursorCharge = precursorCharge;
            this.precursorMass = precursorMass;
            this.labelling = labelling;
            this.peptide = peptide;
            this.theoMass = theoMass;
            this.isDecoy = isDecoy;
            this.score = score;
            this.deltaLCn = deltaLCn;
            this.deltaCn = deltaCn;
            this.matchedPeakNum = matchedPeakNum;
            this.ionFrac = ionFrac;
            this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
            this.explainedAaFrac = explainedAaFrac;
            this.otherPtmPatterns = otherPtmPatterns;
            this.aScore = aScore;
            this.candidates = candidates;
            this.peptideSet = peptideSet;
        }
    }
}
