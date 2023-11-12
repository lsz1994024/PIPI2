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
import proteomics.FM.FMIndex;
import proteomics.FM.FMRes;
import proteomics.Index.BuildIndex;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.ExpAa;
import proteomics.Types.ExpTag;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

import static proteomics.PIPI.minTagLenToReduceProtDb;

public class GetLongTag implements Callable<GetLongTag.Entry> {
    private static final int candisNum = 20;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
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

    public GetLongTag(int scanNum, BuildIndex buildIndex, MassTool massTool, JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock
            , String scanName, int precursorCharge, double precursorMass, SpecProcessor specProcessor, String truth, double ms2Tolerance) {

        this.buildIndex = buildIndex;
        this.massTool = massTool;
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
        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN, ms2Tolerance);

        if (plMap.isEmpty()) {
            return null;
        }

        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);
        Map<String, String> protIdSeqMap = buildIndex.protSeqMap;
        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, minTagLenToReduceProtDb, 99);

        if (!allLongTagList.isEmpty()) {
            Entry entry = new Entry();
            double scoreCutOff = allLongTagList.get(0).getTotalIntensity() - 3;
            int rankCutOff = 10;
            int validTagNum = 0;
            for (int ii = 0; ii < Math.min(rankCutOff, allLongTagList.size()); ii++) {
                if (allLongTagList.get(ii).getTotalIntensity() > scoreCutOff) {
                    validTagNum++;
                } else {
                    break;
                }
            }

            FMIndex fmIndex = buildIndex.fmIndexFull;

            Set<String> usedLongTags = new HashSet<>();
            Queue<ExpTag> tagQueue = new LinkedList<>();

            for (ExpTag longTag : allLongTagList.subList(0, Math.min(validTagNum, allLongTagList.size()))) {
                if (longTag.getTotalIntensity() < longTag.size() * 0.8 + 0.8 && longTag.size() < 10) continue;
                tagQueue.offer(longTag);
            }
            Set<ExpTag> tagThisRound = new HashSet<>();
            Map<String, List<TagRes>> prot_TagResList_Map = new HashMap<>();
            double totalMass = precursorMass + 2*MassTool.PROTON;
            while (!tagQueue.isEmpty()) {
                int ptnForwardCount = 0;
                int ptnBackwardCount = 0;
                FMRes searchForward = null;
                FMRes searchBackward = null;

                ExpTag expTag = tagQueue.poll();
                String tagStr = expTag.getFreeAaString();

                if (!usedLongTags.contains(tagStr)) {
                    usedLongTags.add(tagStr);
                    char[] tagChar = tagStr.toCharArray();
                    searchForward = fmIndex.fmSearch(tagChar);
                    if (searchForward != null) {
                        ptnForwardCount = searchForward.ep - searchForward.sp + 1;
                    }
                }

                ExpTag revExpTag = expTag.revTag(totalMass);
                String revTagStr = revExpTag.getFreeAaString();
                if (!usedLongTags.contains(revTagStr)) {
                    usedLongTags.add(revTagStr);
                    char[] revTagChar = revTagStr.toCharArray();
                    searchBackward = fmIndex.fmSearch(revTagChar);
                    if (searchBackward != null) {
                        ptnBackwardCount = searchBackward.ep - searchBackward.sp + 1;
                    }
                }

                if (ptnForwardCount + ptnBackwardCount > 0) {
                    Set<String> thisRoundForwardProts = new HashSet<>();
                    if (ptnForwardCount > 0) {
                        for (int ii = searchForward.sp; ii <= searchForward.ep; ii++) {
                            int res = Arrays.binarySearch(buildIndex.dotPosArrFull, fmIndex.SA[ii]);
                            thisRoundForwardProts.add(buildIndex.posProtMapFull.get(-res - 2));
                        }
                    }
                    Set<String> ptnForwardGroup = new HashSet<>();
                    for (String prot1 : thisRoundForwardProts) {
                        String seq1 = protIdSeqMap.get(prot1);
                        Set<String> tempGroup = new HashSet<>();
                        tempGroup.add(prot1);
                        for (String prot2 : thisRoundForwardProts) {
                            if (prot1.contentEquals(prot2)) continue;
                            String seq2 = protIdSeqMap.get(prot2);
                            if (seq1.length() != seq2.length()) continue;
                            int errorNum = 0;
                            for (int iChar = 0; iChar < seq1.length(); iChar++) {
                                if (seq1.charAt(iChar) != seq2.charAt(iChar)) {
                                    errorNum++;
                                }
                            }
                            if (errorNum < 0.1 * seq1.length()) {
                                tempGroup.add(prot2);
                            }
                        }
                        if (tempGroup.size() > 1) {
                            ptnForwardGroup.addAll(tempGroup);
                            break;
                        }
                    }

                    Set<String> thisRoundBackwardProts = new HashSet<>();
                    if (ptnBackwardCount > 0) {
                        for (int ii = searchBackward.sp; ii <= searchBackward.ep; ii++) {
                            int res = Arrays.binarySearch(buildIndex.dotPosArrFull, fmIndex.SA[ii]);
                            thisRoundBackwardProts.add(buildIndex.posProtMapFull.get(-res - 2));
                        }
                    }

                    // this is for protein group count, in order to cluster those very alike proteins and make the |P(t)| in paper reasonable
                    Set<String> ptnBackGroup = new HashSet<>();
                    for (String prot1 : thisRoundBackwardProts) {
                        String seq1 = protIdSeqMap.get(prot1);
                        Set<String> tempGroup = new HashSet<>();
                        tempGroup.add(prot1);
                        for (String prot2 : thisRoundBackwardProts) {
                            if (prot1.contentEquals(prot2)) continue;
                            String seq2 = protIdSeqMap.get(prot2);
                            if (seq1.length() != seq2.length()) continue;
                            int errorNum = 0;
                            for (int iChar = 0; iChar < seq1.length(); iChar++) {
                                if (seq1.charAt(iChar) != seq2.charAt(iChar)) {
                                    errorNum++;
                                }
                            }
                            if (errorNum < 0.1 * seq1.length()) {
                                tempGroup.add(prot2);
                            }
                        }
                        if (tempGroup.size() > 1) {
                            ptnBackGroup.addAll(tempGroup);
                            break;
                        }
                    }

                    if (ptnForwardCount > 0) {
                        List<Double> normedIaaList = new ArrayList<>(expTag.expAaList.size());
                        for (ExpAa aa : expTag.expAaList) {
                            normedIaaList.add(aa.getTotalIntensity()*0.5/(ptnForwardCount - ptnForwardGroup.size() + 1));
                        }
                        for (int ii = searchForward.sp; ii <= searchForward.ep; ii++) {
                            int absPos = fmIndex.SA[ii];
                            int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrFull, absPos);
                            String prot = buildIndex.posProtMapFull.get(dotIndex);
                            int relPos = absPos - buildIndex.dotPosArrFull[dotIndex] - 1;
                            if (prot_TagResList_Map.containsKey(prot)) {
                                prot_TagResList_Map.get(prot).add(new TagRes(tagStr, relPos, normedIaaList));
                            } else {
                                List<TagRes> tmpTagResList = new LinkedList<>();
                                tmpTagResList.add(new TagRes(tagStr, relPos, normedIaaList));
                                prot_TagResList_Map.put(prot, tmpTagResList);
                            }
                        }
                    }
                    if (ptnBackwardCount > 0) {
                        List<Double> normedIaaList = new ArrayList<>(revExpTag.expAaList.size());
                        for (ExpAa aa : revExpTag.expAaList) {
                            normedIaaList.add(aa.getTotalIntensity()*0.5/(ptnBackwardCount - ptnBackGroup.size() + 1));
                        }
                        for (int ii = searchBackward.sp; ii <= searchBackward.ep; ii++) {
                            int absPos = fmIndex.SA[ii];
                            int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrFull, absPos);
                            String prot = buildIndex.posProtMapFull.get(dotIndex);
                            int relPos = absPos - buildIndex.dotPosArrFull[dotIndex] - 1;
                            if (prot_TagResList_Map.containsKey(prot)) {
                                prot_TagResList_Map.get(prot).add(new TagRes(revTagStr, relPos, normedIaaList));
                            } else {
                                List<TagRes> tmpTagResList = new LinkedList<>();
                                tmpTagResList.add(new TagRes(revTagStr, relPos, normedIaaList));
                                prot_TagResList_Map.put(prot, tmpTagResList);
                            }
                        }
                    }
                    break; // if it is processing original tags, dont break
                }
                if (tagStr.length() > 6) tagThisRound.add(expTag);

//                if (tagQueue.isEmpty()) {
//                    List<ExpTag> tempList = new LinkedList<>();
//                    for (ExpTag tag : tagThisRound) {
//                        ExpTag tagL = tag.subTag(0, tag.size() - 1);
//                        if (!usedLongTags.contains(tagL.getFreeAaString())) {
//                            tempList.add(tagL);
//                        }
//                        ExpTag tagR = tag.subTag(1, tag.size());
//                        if (!usedLongTags.contains(tagR)) {
//                            tempList.add(tagR);
//                        }
//                    }
//                    tempList.sort(Comparator.comparingDouble(ExpTag::getTotalIntensity));
//                    for (int j = tempList.size() - 1; j >= 0; j--) {
//                        tagQueue.offer(tempList.get(j));
//                    }
//                }
            }

//            System.out.println();
            entry.prot_TagResList_Map = prot_TagResList_Map;
            entry.scanName = this.scanName;
            return entry;
        } else {
            return null;
        }
    }

    public class Entry {
        public Map<String, List<TagRes>> prot_TagResList_Map = new HashMap<>();
        public String scanName;
        Entry() {
        }
    }
    public class TagRes {
        public String tagSeq;
        public int relPos;
        public List<Double> normedIaaList = null;
        public TagRes(String tagSeq, int relPos, List<Double> normedIaaList) {
            this.relPos = relPos;
            this.tagSeq = tagSeq;
            this.normedIaaList = normedIaaList;
        }
    }
}
