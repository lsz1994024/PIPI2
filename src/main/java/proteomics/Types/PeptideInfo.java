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
import java.util.HashSet;
import java.util.Set;

public class PeptideInfo implements Cloneable{
    public String freeSeq;
    public boolean isDecoy;
    public Set<String> protIdSet = new HashSet<>();
    public char leftFlank;
    public char rightFlank;


    public PeptideInfo(String freeSeq, boolean isDecoy, char leftFlank, char rightFlank) {
        this.isDecoy = isDecoy;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
        this.freeSeq = freeSeq;
    }

    public PeptideInfo clone() throws CloneNotSupportedException {
        super.clone();
        PeptideInfo other = new PeptideInfo(freeSeq, isDecoy, leftFlank, rightFlank);
        other.protIdSet.addAll(this.protIdSet);
        return other;
    }

}
