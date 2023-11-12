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
import proteomics.Index.BuildIndex;

import java.util.*;
import java.util.concurrent.Callable;

public class BuildDecoyProts implements Callable<BuildDecoyProts.Entry> {
    private Map<String, String> parameterMap;
//    private Set<String> reducedProtIdSet;

    private BuildIndex buildIndex;

    private String protId;
    public BuildDecoyProts(Map<String, String> parameterMap, BuildIndex buildIndex, String protId) throws Exception {
        this.parameterMap = parameterMap;
        this.buildIndex = buildIndex;
        this.protId = protId;
    }

    @Override
    public Entry call() throws Exception {
        boolean addDecoy = parameterMap.get("add_decoy").contentEquals("1");

        Entry entry = new Entry();
        entry.protId = protId;
        String protSeq = buildIndex.protSeqMap.get(protId).replace('I', 'L');
        if (addDecoy) {
            String decoyProtSeq = DbTool.shuffleProtKeepKR(protSeq, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), Integer.valueOf(parameterMap.get("is_from_C_term_1")) == 1).replace('I', 'L');
            entry.decoyProtSeq = decoyProtSeq;
        }
        return entry;
    }

    public class Entry {
        String protId = "";
        String decoyProtSeq = "";
        Entry() {
        }
    }
}
