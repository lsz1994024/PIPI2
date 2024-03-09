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

package proteomics.Segment;


import ProteomicsLibrary.DbTool;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Types.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static proteomics.PIPI.MIN_PEAK_SUM_INFER_AA;

public class InferSegment {
    private static final Logger logger = LoggerFactory.getLogger(InferSegment.class);
    private static final int minTagNum = 200;
    private static final int regionNum = 10;
    private static final int topNumInEachRegion = 20;
    public static final byte N_TAG = -1;
    public static final byte NON_NC_TAG = 0;
    public static final byte C_TAG = 1;
    private static final Pattern pattern = Pattern.compile("([nc][0-9a-i])?([A-Z#$].?)");
    private final double ms2Tolerance;
    private TreeMap<Segment, Integer> aaVectorTemplate = new TreeMap<>();
    private Map<Double, String> augedMassAaMap = new HashMap<>(35, 1);// augmented amino acid map. Normal aa plus aa with mod
    private double maxAugedMass = 0;
    private final Double[] augedMassArray;
    private Map<String, Double> extraAaMassMap = new HashMap<>(35, 1);
    private double[] nTermPossibleMod = null;
    private double[] cTermPossibleMod = null;
    private MassTool massTool;
    public String aaModLabel = "~!@#$%^&*";
    public String nModLabel = "qwertyuio";
    public String cModLabel = "asdfghjkl";

    public InferSegment(MassTool massTool, Map<String, String> parameterMap, Map<Character, Double> fixModMap) throws Exception {
        this.massTool = massTool;
        this.ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        Map<Character, Double> massTable = massTool.getMassTable();

        char[] standardAaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W'};

        Map<Double, Character> oriMassAaMap = new HashMap<>(25, 1);
        for (char aa : standardAaArray) {
            // # = I/L.
            oriMassAaMap.put(massTable.get(aa), aa);
        }

        Character[] aaArray = oriMassAaMap.values().toArray(new Character[0]);

        for (char aa1 : aaArray) {
            for (char aa2 : aaArray) {
                for (char aa3 : aaArray) {
                    aaVectorTemplate.put(new Segment(String.format(Locale.US, "%c%c%c", aa1, aa2, aa3)), 0);
                }
            }
        }

        int idx = 0;
        for (Segment segment : aaVectorTemplate.keySet()) {
            aaVectorTemplate.put(segment, idx);
            ++idx;
        }

        // generate a mass aa map containing modified amino acid
        for (double k : oriMassAaMap.keySet()) {
            augedMassAaMap.put(k, oriMassAaMap.get(k).toString());
        }
        int iLabelAa = 0;
        int iLabelN = 0;
        int iLabelC = 0;
        for (String k : parameterMap.keySet()) {
            if (!k.startsWith("mod")) continue;

            String v = parameterMap.get(k);
            if (v.startsWith("0.0")) break;

            String[] modStr = v.split(",");
            double modMass = Double.valueOf(modStr[0]);
            char modSite = modStr[1].charAt(0);
            char label;
            int modPosition = Integer.valueOf(modStr[2]);
            VarPtm varPtmByUser = new VarPtm(modMass, modSite, modPosition, modStr[3], "ByUser", 1);
            if (modPosition == 4) {//  position anywhere, highest prority
                boolean added = false;
                label = aaModLabel.charAt(iLabelAa);
                double augedMass = massTable.get(modSite) + modMass;
                boolean shouldAdd = true;
                for (double tempMass : augedMassAaMap.keySet()) {
                    if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                        logger.warn(String.format(Locale.US, "%s is bad. %s and %s have conflict mass values(%f vs %f).", v, modSite, augedMassAaMap.get(tempMass), augedMass, tempMass));
                        shouldAdd = false;
                    }
                }
                if (Math.abs(fixModMap.get(modSite)) < 0.1 && shouldAdd) {
                    String augedAa = ""+modSite+label;
                    augedMassAaMap.put(augedMass, augedAa);
                    extraAaMassMap.put(augedAa, modMass);
                    added = true;
//                    massTool.labelMassMap.put(label, modMass);
                }

                if (added) {
                    iLabelAa++;
                    massTool.labelVarPtmMap.put(label, varPtmByUser);
                }

            } else if (modPosition == 0 || modPosition == 2) {
                label = nModLabel.charAt(iLabelN);
                boolean added = false;
                if (modSite == 'X') { //on any aa
                    for (char oriAa : aaArray){
                        double augedMass = massTable.get(oriAa) + modMass;
                        boolean shouldAdd = true;
                        for (double tempMass : augedMassAaMap.keySet()) {
                            if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                                logger.warn(String.format(Locale.US, "%s is bad. %s and %s have conflict mass values(%f vs %f).", v, oriAa, augedMassAaMap.get(tempMass), augedMass, tempMass));
                                shouldAdd = false;
                            }
                        }
                        if (Math.abs(fixModMap.get(oriAa)) < 0.1 && shouldAdd) {
                            String augedAa = ""+oriAa+label;
                            augedMassAaMap.put(augedMass, augedAa);
                            extraAaMassMap.put(augedAa, modMass);
                            added = true;
                        }
                    }
                } else { //on single aa
                    double augedMass = massTable.get(modSite) + modMass;
                    boolean shouldAdd = true;
                    for (double tempMass : augedMassAaMap.keySet()) {
                        if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                            logger.warn(String.format(Locale.US, "%s and %s have conflict mass values(%f vs %f).", v, augedMassAaMap.get(tempMass), augedMass, tempMass));
                            shouldAdd = false;
                        }
                    }
                    if (Math.abs(fixModMap.get(modSite)) < 0.1 && shouldAdd) {
                        String augedAa = ""+modSite+label;
                        augedMassAaMap.put(augedMass, augedAa);
                        extraAaMassMap.put(augedAa, modMass);
                        added = true;
                    }
                }

                if (added) {
                    iLabelN++;
                    massTool.labelVarPtmMap.put(label, varPtmByUser);
                }
            } else {// position C term, lowest  prority
                label = cModLabel.charAt(iLabelC);
                boolean added = false;
                if (modSite == 'X') { //on any aa
                    for (char oriAa : aaArray){
                        double augedMass = massTable.get(oriAa) + modMass;
                        boolean shouldAdd = true;
                        for (double tempMass : augedMassAaMap.keySet()) {
                            if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                                logger.warn(String.format(Locale.US, "%s is bad. %s and %s have conflict mass values(%f vs %f).", v, oriAa, augedMassAaMap.get(tempMass), augedMass, tempMass));
                                shouldAdd = false;
                            }
                        }
                        if (Math.abs(fixModMap.get(oriAa)) < 0.1 && shouldAdd) {
                            String augedAa = ""+oriAa+label;
                            augedMassAaMap.put(augedMass, augedAa);
                            extraAaMassMap.put(augedAa, modMass);
                            added = true;
                        }
                    }
                } else {
                    double augedMass = massTable.get(modSite) + modMass;
                    boolean shouldAdd = true;
                    for (double tempMass : augedMassAaMap.keySet()) {
                        if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                            logger.warn(String.format(Locale.US, "%s and %s have conflict mass values(%f vs %f).", v, augedMassAaMap.get(tempMass), augedMass, tempMass));
                            shouldAdd = false;
                        }
                    }
                    if (Math.abs(fixModMap.get(modSite)) < 0.1 && shouldAdd) {
                        String augedAa = ""+modSite+label;
                        augedMassAaMap.put(augedMass, augedAa);
                        extraAaMassMap.put(augedAa, modMass);
                        added = true;
                    }
                }
                if (added) {
                    iLabelC++;
                    massTool.labelVarPtmMap.put(label, varPtmByUser);
                }
            }
        }
        augedMassArray = augedMassAaMap.keySet().toArray(new Double[0]);
        maxAugedMass = Collections.max(augedMassAaMap.keySet()) + ms2Tolerance;
    }

    public SparseBooleanVector generateSegmentBooleanVector(String peptide) {

        String normalizedPeptide = normalizeSequence(DbTool.getSequenceOnly(peptide));

        Set<Integer> tempSet = new HashSet<>(DbTool.getSequenceOnly(peptide).length() + 1, 1);
        for (int i = 0; i <= normalizedPeptide.length() - 3; ++i) {
            tempSet.add(aaVectorTemplate.get(new Segment(normalizedPeptide.substring(i, i + 3))));
        }
        return new SparseBooleanVector(tempSet);
    }
    public static String normalizeSequence(String seq) {
        return seq.replaceAll("I", "L");
    }

    public List<ExpTag> cleanAbundantTagsPrefix(List<ExpTag> allLongTagsList, int minTagLen) {
        if (allLongTagsList.isEmpty()) return allLongTagsList;
        List<ExpTag> cleanTagList = new LinkedList<>();
        Set<String> prefixSet = new HashSet<>();
        for (ExpTag tag : allLongTagsList) {
            String prefix = tag.getFreeAaString().substring(0,minTagLen);
            boolean shouldAdd = true;
            if (!prefixSet.contains(prefix)) {
                prefixSet.add(prefix);
            } else {
                shouldAdd = false;
            }
            if (shouldAdd) {
                cleanTagList.add(tag);
            }
        }
        return cleanTagList;
    }


    private static int minDistance(String word1, String word2) {

        int dp[][] = new int[word1.length() + 1][word2.length() + 1];

        for (int i = 0; i < word1.length() + 1; i++) {
            dp[i][0] = i;
        }
        for (int i = 0; i < word2.length() + 1; i++) {
            dp[0][i] = i;
        }

        for (int i = 1; i < word1.length() + 1; i++) {
            for (int j = 1; j < word2.length() + 1; j++) {
                if (word1.charAt(i - 1) == word2.charAt(j - 1)) {
                    dp[i][j] = dp[i - 1][j - 1];
                } else {
                    dp[i][j] = 1 + Math.min(dp[i - 1][j - 1], Math.min(dp[i - 1][j], dp[i][j - 1]));
                }
            }
        }
        return dp[word1.length()][word2.length()];
    }
    public List<ExpTag> getLongTag(TreeMap<Double, Double> plMap, double cTermMz, int scanNum, int minTagLenToExtractLocal, int maxTagLenToExtractLocal) throws Exception {
        Double[] mzArray = plMap.keySet().toArray(new Double[0]);
        Double[] intensityArray = plMap.values().toArray(new Double[0]);
        List<ExpTag> outputList = new LinkedList<>();

        Set<Pair<Integer, Integer>> edgeSet = new HashSet<>();
        Set<Integer> nodeSet = new HashSet<>();
        Set<Integer> startNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
        Set<Integer> endNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
        Map<Pair<Integer, Integer>, ExpAa> edgeInfoMap = new HashMap<>();
        for (int i = 0; i < mzArray.length - 1; ++i) {
            double mz1 = mzArray[i];
            int isNorC = NON_NC_TAG;
            if (Math.abs(mz1 - MassTool.PROTON) <= ms2Tolerance) {
                isNorC = N_TAG;
            } else if (Math.abs(mz1 - MassTool.PROTON - massTool.H2O) <= ms2Tolerance) {
                isNorC = C_TAG;
            }
            double intensity1 = intensityArray[i];
            if (intensity1 < 0.1) continue;//todo

            for (int j = i + 1; j < mzArray.length; ++j) {
                double mz2 = mzArray[j];
                if (mz2 > mz1 + maxAugedMass + 5) break; // no need to try further because the can not match to any aa
                double intensity2 = intensityArray[j];
                if (intensity2 < 0.1) continue;//todo

                if ( (intensity1 + intensity2) < MIN_PEAK_SUM_INFER_AA) continue;  //todo

                String aa = inferAA(mz1, mz2, isNorC);

                if (aa != null ) {
                    Pair<Integer, Integer> edge = new Pair<>(i,j);
                    edgeSet.add(edge);
                    nodeSet.add(i);
                    nodeSet.add(j);
                    startNodeSet.remove(j);
                    endNodeSet.remove(i);
                    edgeInfoMap.put(edge, new ExpAa(aa, aa.charAt(0), mz1, mz2, intensity1, intensity2, isNorC));
                }
            }
        }
        if (edgeInfoMap.size() > 120){ //only use the top 100 edges
            List<Map.Entry<Pair<Integer, Integer>, ExpAa>> edgeInfoList = new ArrayList<>(edgeInfoMap.entrySet());
            Collections.sort(edgeInfoList, Comparator.comparing(o -> o.getValue().getTotalIntensity(), Comparator.reverseOrder()));
            nodeSet.clear();
            edgeSet.clear();
            startNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
            endNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
            int i = 0;
            for (Map.Entry<Pair<Integer, Integer>, ExpAa> entry : edgeInfoList){
                if (i > 120) {
                    edgeInfoMap.remove(entry.getKey());
                } else {
                    edgeSet.add(entry.getKey());
                    nodeSet.add(entry.getKey().getFirst());
                    nodeSet.add(entry.getKey().getSecond());
                    startNodeSet.remove(entry.getKey().getSecond());
                    endNodeSet.remove(entry.getKey().getFirst());
                }
                i++;
            }

        }

        startNodeSet.retainAll(nodeSet);
        endNodeSet.retainAll(nodeSet);
        Graph g = new Graph(edgeSet, nodeSet);
        ArrayList<ArrayList<Integer>> allPath = g.getAllPaths(startNodeSet, endNodeSet);
        for (ArrayList<Integer> path : allPath) {
            if (path.size() < minTagLenToExtractLocal+1) continue; //aa length is 6 . peaks number is 7.
            if (path.size() > maxTagLenToExtractLocal + 1) {
                for (int i = 0; i < path.size()-maxTagLenToExtractLocal; i++) {
                    List<ExpAa> expAaList = new ArrayList<>();
                    for (int ii = i; ii < i+maxTagLenToExtractLocal; ii++){
                        int j = ii + 1;
                        expAaList.add(edgeInfoMap.get(new Pair<>(path.get(ii),path.get(j))));
                    }
                    outputList.add(new ExpTag(expAaList));
                }
            } else {
                List<ExpAa> expAaList = new ArrayList<>();
                for (int i = 0; i < path.size()-1; i++){
                    int j = i + 1;
                    expAaList.add(edgeInfoMap.get(new Pair<>(path.get(i),path.get(j))));
                }
                outputList.add(new ExpTag(expAaList));
            }
        }
        outputList.sort(Comparator.comparingDouble(ExpTag::getTotalIntensity).reversed());
        boolean[] shouldKeep = new boolean[outputList.size()];
        Set<String> addedTags = new HashSet<>();
        int i = 0;
        for(ExpTag tag : outputList) {
            if (addedTags.contains(tag.getFreeAaString())) {
                shouldKeep[i] = false;
            } else {
                shouldKeep[i] = true;
                addedTags.add(tag.getFreeAaString());
            }
            i++;
        }
        List<ExpTag> finalList = new LinkedList<>();
        for (int j = 0; j < outputList.size(); j++){
            if (shouldKeep[j]) {
                finalList.add(outputList.get(j));
            }
        }
        return finalList;
    }

    public int getLongTagNumForProt(String prot) {
        String normalizedProt = normalizeSequence(DbTool.getSequenceOnly(prot));
        Set<String> tagSet = new HashSet<>();
        for (int i = 0; i <= normalizedProt.length() - 7; ++i) {
            Segment seg = new Segment(normalizedProt.substring(i, i + 7));
            String tag = seg.toString();
            tagSet.add(tag);
        }

        return tagSet.size();
    }

    private String inferAA(double mz1, double mz2, int isNorC) {
        double mzDiff = mz2 - mz1;
        String aa = null;
        for (Map.Entry<Double, String> massAa : augedMassAaMap.entrySet()) {
            String augedAa = massAa.getValue();
            if ((isNorC == N_TAG && cModLabel.contains(augedAa.substring(augedAa.length()-1)))  // use n peak then should no C mod
                    || (isNorC == C_TAG && nModLabel.contains(augedAa.substring(augedAa.length()-1))) // use c peak then should no N mod
                    || (isNorC == NON_NC_TAG && cModLabel.contains(augedAa.substring(augedAa.length()-1))) // use non-nc peak then should no C mod
                    || (isNorC == NON_NC_TAG && nModLabel.contains(augedAa.substring(augedAa.length()-1))) // use non-nc peak then should no N mod
            ) {
                continue;
            }

            if (Math.abs(mzDiff - massAa.getKey()) <= 1 * ms2Tolerance) {
                aa = massAa.getValue();// including M~
            }
        }

        if (aa != null && isNorC == C_TAG && (!"KR".contains(aa))) {
            aa = null;
        }
        return aa;
    }

    public TreeMap<Double, Double> addVirtualPeaks(double precursorMass, TreeMap<Double, Double> plMap) {
        double totalMass = precursorMass + 2 * MassTool.PROTON;
        TreeMap<Double, Double> finalPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            finalPlMap.put(mz, plMap.get(mz));
        }
        for (double mz : plMap.keySet()) {
            double anotherMz = totalMass - mz;
            double leftMz = anotherMz - ms2Tolerance;
            double rightMz = anotherMz + ms2Tolerance;
            NavigableMap<Double, Double> temp = null;
            temp = plMap.subMap(leftMz, true, rightMz, true);

        }

        finalPlMap.put(MassTool.PROTON, 1d);
        double cTermMz = precursorMass - massTool.H2O + MassTool.PROTON;
        double leftMz = cTermMz - ms2Tolerance;
        double rightMz = cTermMz + ms2Tolerance;
        NavigableMap<Double, Double> temp = null;
        try {
            temp = plMap.subMap(leftMz, true, rightMz, true);
        } catch (IllegalArgumentException ex) {}
        if ((temp == null) || (temp.isEmpty())) {
            finalPlMap.put(cTermMz, 1d);
        }
        finalPlMap.put(MassTool.PROTON + massTool.H2O, 1d);

        return finalPlMap;
    }
}
