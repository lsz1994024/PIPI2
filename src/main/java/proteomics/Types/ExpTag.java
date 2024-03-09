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

package proteomics.Types;



import java.util.ArrayList;
import java.util.List;

public class ExpTag implements Comparable<ExpTag> {
    public int isNorC = 0;
    private final ExpAa[] expAaArray;
    private int hashCode;
    private final double totalIntensity;
    private final String freeAaString;
    private final String ptmAaString;
    private int regionIdx;
    public List<ExpAa> expAaList = null;
    public double[] intensityArray;
    public ExpTag(ExpAa aa1, ExpAa aa2, ExpAa aa3) {
        expAaArray = new ExpAa[]{aa1, aa2, aa3};
        String toString = expAaArray[0].toString() + "-" + expAaArray[1].toString() + "-" + expAaArray[2].toString();
        hashCode = toString.hashCode();

        StringBuilder sbFreeSeq = new StringBuilder(expAaArray.length);
        StringBuilder sbPtmSeq = new StringBuilder(2* expAaArray.length - 1);
        for (ExpAa aa : expAaArray) {
            sbFreeSeq.append(aa.getPtmFreeAA());
            sbPtmSeq.append(aa.getAA());
        }
        freeAaString = sbFreeSeq.toString();
        ptmAaString = sbPtmSeq.toString();

        double intensity = expAaArray[0].getHeadIntensity();
        for (ExpAa aa : expAaArray) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }
    public ExpTag(ExpAa aa1, ExpAa aa2, ExpAa aa3, ExpAa aa4) {
        expAaArray = new ExpAa[]{aa1, aa2, aa3, aa4};
        String toString = expAaArray[0].toString() + "-" + expAaArray[1].toString() + "-" + expAaArray[2].toString()+ "-" + expAaArray[3].toString();
        hashCode = toString.hashCode();

        StringBuilder sbFreeSeq = new StringBuilder(expAaArray.length);
        StringBuilder sbPtmSeq = new StringBuilder(2* expAaArray.length - 1);
        for (ExpAa aa : expAaArray) {
            sbFreeSeq.append(aa.getPtmFreeAA());
            sbPtmSeq.append(aa.getAA());
        }
        freeAaString = sbFreeSeq.toString();
        ptmAaString = sbPtmSeq.toString();

        double intensity = expAaArray[0].getHeadIntensity();
        for (ExpAa aa : expAaArray) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }

    public ExpTag(List<ExpAa> expAaList) {
        isNorC = expAaList.get(0).isNorC;
        intensityArray = new double[expAaList.size()+1];
        this.expAaList = expAaList;
        intensityArray[0] = expAaList.get(0).getHeadIntensity();
        for (int i = 0; i < expAaList.size(); i++) {
            intensityArray[i+1] = expAaList.get(i).getTailIntensity();
        }
        expAaArray = expAaList.toArray(new ExpAa[0]);
        String toString = expAaArray[0].toString();
        for (int i = 1; i < expAaArray.length; i++){
            toString += "-" + expAaArray[i].toString();
        }
        hashCode = toString.hashCode();
        StringBuilder sbFreeSeq = new StringBuilder(expAaArray.length);
        StringBuilder sbPtmSeq = new StringBuilder(2* expAaArray.length - 1);
        for (ExpAa aa : expAaArray) {
            sbFreeSeq.append(aa.getPtmFreeAA());
            sbPtmSeq.append(aa.getAA());
        }
        freeAaString = sbFreeSeq.toString();
        ptmAaString = sbPtmSeq.toString();

        double intensity = expAaArray[0].getHeadIntensity();
        for (ExpAa aa : expAaArray) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        return (other instanceof ExpTag) && (this.hashCode() == other.hashCode());
    }

    public boolean approximateEquals(ExpTag other, double tolerance) {
        for (int i = 0; i < this.size(); ++i) {
            if (!this.get(i).approximateEquals(other.get(i), tolerance)) {
                return false;
            }
        }
        return true;
    }

    public int compareTo(ExpTag other) {
        return Double.compare(expAaArray[0].getHeadLocation(), other.expAaArray[0].getHeadLocation());
    }

    public String getFreeAaString() {
        return freeAaString;
    }
    public String getPtmAaString() {
        return ptmAaString;
    }

    public String toString() {
        return isNorC+","+ ptmAaString ;
    }

    public double getTotalIntensity() {
        return totalIntensity;
    }

    public double getHeadLocation() {
        return expAaArray[0].getHeadLocation();
    }

    public double getTailLocation() {
        return expAaArray[expAaArray.length - 1].getTailLocation();
    }

    public ExpTag clone() throws CloneNotSupportedException {
        super.clone();
        if (expAaArray.length == 3 ){
            return new ExpTag(expAaArray[0].clone(), expAaArray[1].clone(), expAaArray[2].clone());
        }else {
            return new ExpTag(expAaArray[0].clone(), expAaArray[1].clone(), expAaArray[2].clone(), expAaArray[3].clone());
        }
    }

    public int size() {
        return expAaArray.length;
    }

    public ExpAa get(int i) {
        return expAaArray[i];
    }


    public ExpTag revTag(double totalMass){
        List<ExpAa> revExpAaList = new ArrayList<>(expAaArray.length);
        for (int i = expAaArray.length-1; i>=0; i--) {
            revExpAaList.add(expAaArray[i].revAA(totalMass));
        }
        ExpTag revedTag = new ExpTag(revExpAaList);
        revedTag.isNorC = isNorC;
        return revedTag;
    }

    public ExpTag subTag(int start, int end){ // start included  end not included
        List<ExpAa> subExpAaList = new ArrayList<>(end-start);
        for (int i = start; i < end; i++) {
            subExpAaList.add(expAaArray[i]);
        }
        ExpTag subTag = new ExpTag(subExpAaList);
        subTag.isNorC = isNorC;
        return subTag;
    }

}
