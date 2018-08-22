import java.awt.List;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class NeedlemanWunsch {
    
    public static void main (String[] args) throws IOException {
    	
        String  FILENAME = args[0];
    	ArrayList lines =  (ArrayList) Files.readAllLines(Paths.get(FILENAME));
    	HashMap<ArrayList[], Integer> map = new HashMap<ArrayList[], Integer>();
    	
    	// get the dp tables
        int [][] dpTable1 = createDPTable((String)lines.get(0),(String)lines.get(1)); // 1,2
        int [][] dpTable2 = createDPTable((String)lines.get(0),(String)lines.get(2)); // 1,3
        int [][] dpTable3 = createDPTable((String)lines.get(1),(String)lines.get(2)); // 2,3
        
        // get the scores of each alignment
        int score1 = getScore(dpTable1);
        int score2 = getScore(dpTable2);
        int score3 = getScore(dpTable3);
        
        // do the alignment for each pair
        ArrayList[] list1 = NWPairwiseAlignment(dpTable1,(String)lines.get(0),(String)lines.get(1)); // 1,2
        ArrayList[] list2 = NWPairwiseAlignment(dpTable2,(String)lines.get(0),(String)lines.get(2)); // 1,3
        ArrayList[] list3 = NWPairwiseAlignment(dpTable3,(String)lines.get(1),(String)lines.get(2)); // 2,3
        
        // print alignment and corresponding scores
        System.out.println("Alignment 1&2 score: " + score1);
        printSeq(list1[0]);
        printSeq(list1[1]);
        System.out.println();
        System.out.println("Alignment 1&3 score: " + score2);
        printSeq(list2[0]);
        printSeq(list2[1]);
        System.out.println();
        System.out.println("Alignment 2&3 score: " + score3);
        printSeq(list3[0]);
        printSeq(list3[1]);
        System.out.println();
        

        // store the pairs and their corresponding scores
        map.put(list1, score1);
        map.put(list2, score2);
        map.put(list3, score3);
        
        // get the alignment with the maximum score.
        Map.Entry<ArrayList[], Integer> maxEntry = null;
        for (Map.Entry<ArrayList[], Integer> entry : map.entrySet()) {
          if (maxEntry == null || entry.getValue() > maxEntry.getValue()) {
            maxEntry = entry;
          }
        }
        ArrayList[] maxList = maxEntry.getKey(); 
        
        String intermediateSeq = mergeTwoSeq(maxList[0],maxList[1]);
        
        int [][] dpTableMax = createDPTable((String)lines.get(1),intermediateSeq); 
        int scoreMax = getScore(dpTableMax);
        ArrayList[] resultList = NWPairwiseAlignment(dpTableMax,(String)lines.get(1),intermediateSeq);
        
        // show triple alignment 
        String seq1 = convert(maxList[0]);
        String seq2 = convert(maxList[1]);
        String seq3 = convert(resultList[0]);
        seq1 = handleResidual(seq1,intermediateSeq);
        seq2 = handleResidual(seq2,intermediateSeq);
        System.out.println("The maximum alignment score of the final multiple sequence alignment: " + score3 );
        System.out.println("The final multiple alignment achieving this maximum score:");
        System.out.println(seq1);
        System.out.println(seq2);
        System.out.println(seq3);
        
    }
    
    public static int[][] createDPTable(String seq1, String seq2){
        
        // create DP table
        int [][] dpTable = new int[seq1.length()+1][seq2.length()+1];
        
        // initialize the first row
        for (int i = 0; i< dpTable.length; i++)
            dpTable[i][0] = 0; // penalty is 0
        
        // initialize the first column
        for (int i = 0; i< dpTable[0].length; i++)
            dpTable[0][i] = 0; // penalty is 0
        
        for(int i = 1; i<dpTable.length; i++){
            for(int j = 1; j<dpTable[0].length;j++){
                
                int up = dpTable[i-1][j];
                int left = dpTable[i][j-1];
                int diag = dpTable[i-1][j-1];
                
                if (seq1.charAt(i-1) == seq2.charAt(j-1))
                    dpTable[i][j] = Math.max(Math.max(up,left),diag+1);
                else
                    dpTable[i][j] = Math.max(Math.max(up,left),diag);
            }
        }
        
        // score is the last element3 of the table.
        return dpTable;
    }
    
    public static ArrayList[] NWPairwiseAlignment(int[][] dpTable, String seq1, String seq2){
        
        int i = seq1.length();
        int j = seq2.length();
        
        int cnt1 = seq1.length()-1;
        int cnt2 = seq2.length()-1;
        
        // safe enough space
        ArrayList firstSeq = new ArrayList();
        ArrayList secondSeq = new ArrayList();
        
        // start form bottom right
        while(i>0 && j>0){
            
            int left = dpTable[i-1][j];
            int up = dpTable[i][j-1];
            int diag = dpTable[i-1][j-1];
            
            // go to up diagonal
             if ((dpTable[i][j]==diag+1)){
                firstSeq.add(seq1.charAt(cnt1));
                secondSeq.add(seq2.charAt(cnt2));
                i--;
                j--;
                cnt1--;
                cnt2--;
            }
            
            //go to up
            else if(dpTable[i][j]==up){
                firstSeq.add('-'); // empty
                secondSeq.add(seq2.charAt(cnt2));
                j--;
                cnt2--;
            }
            
            // go to left
            else if(dpTable[i][j]==left){
                firstSeq.add(seq1.charAt(cnt1));
                secondSeq.add('-');
                i--;
                cnt1--;
            }
           
        }
        
        // since we AND the index, we need to make sure that we reach top left
        for (;i>0; i--){
            firstSeq.add(seq1.charAt(cnt1)); 
            secondSeq.add('-');
            cnt1--;
        }
        for (;j>0; j--){
            firstSeq.add('-'); // empty
            secondSeq.add(seq2.charAt(cnt2));
            cnt2--;
        }
         
        Collections.reverse(firstSeq);
        Collections.reverse(secondSeq);
        
        ArrayList [] alignedSeqs = new ArrayList[2];
        alignedSeqs[0] = firstSeq;
        alignedSeqs[1] = secondSeq;
        
        return alignedSeqs;
    }
    
    public static int getScore(int[][] dpTable){
    	return dpTable[dpTable.length-1][dpTable[0].length-1];
    }
    
    public static String mergeTwoSeq(ArrayList a1, ArrayList a2){
    
    	String s = "";
    	for(int i = 0; i<a1.size(); i++){
    		if(a1.get(i) == a2.get(i))
    			s=s+a1.get(i);
    		else if((char) a1.get(i) == '-')
    			s=s+a2.get(i);
    		else
    			s=s+a1.get(i);
    	}
    	
    	return s;
    }
    
    public static String convert(ArrayList a1){
        
    	String s = "";
    	for(int i = 0; i<a1.size(); i++){
    			s=s+a1.get(i);	
    	}
    	return s;
    }
    
    
    
    public static void printSeq(ArrayList a1){
    	 for(int i = 0 ; i<a1.size(); i++){
         	System.out.print(a1.get(i));
         }
    	 System.out.println();
    	 
    }
    
    public static String handleResidual(String a, String b){
    	int i =  b.length()-a.length() ;
    	while (i>=0){
    		a = a+'-';
    		i--;
    	}
    	return a;
    }
    
    public static String read(String filename) throws Exception, Exception {
    	try(BufferedReader br = new BufferedReader(new FileReader("file.txt"))) {
    	    StringBuilder sb = new StringBuilder();
    	    String line = br.readLine();

    	    while (line != null) {
    	        sb.append(line);
    	        sb.append(System.lineSeparator());
    	        line = br.readLine();
    	    }
    	    String everything = sb.toString();
    	}
    	return "";
    } 
}
