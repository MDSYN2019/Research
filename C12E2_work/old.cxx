
/*
 
  void AllocationVec() {
      The way that the C12E2 and C12E2-M batches are allocated in the simulation follows
      the fact that this was directly replicated in LAMMPS; hence, there isn't a single 
      block of list, but rather 4 batches of C12E2 and C12E2-M portions that we need to 
      allocate 
      
      Each batch has 250 molecules, which is why we loop from 0 to 249 
    
    for (int i = 0; i <= 249; i++) {

      // C12E2 batches

      int batch1 = 0;
      int batch2 = 3500;
      int batch3 = 7000;
      int batch4 = 10500;

      // C12E2-M batches

      int batch1M = 1750;
      int batch2M = 5250;	
      int batch3M = 8750;
      int batch4M = 12250;
      
      // Allocate index
      

          
      A7 = 7*(i);
      A6_1= 7*(i) + 1;
      A6_2= 7*(i) + 2;
      A3_1 = 7*(i) + 3;
      A3_2 = 7*(i) + 4;
      A3_3 = 7*(i) + 5;
      A4 = 7*(i) + 6;
    
      // Second batch
    
      A7_2 = 7*(i) + 3501;
      A6_1_2 = 7*(i) + 3502;
      A6_2_2 = 7*(i) + 3503;
      A3_1_2 = 7*(i) + 3504;
      A3_2_2 = 7*(i) + 3505;
      A3_3_2 = 7*(i) + 3506;
      A4_2 = 7*(i) + 3507;
    
      // Third batch
      
      A7_3 = 7*(i) + 7001;
      A6_1_3 = 7*(i) + 7002;
      A6_2_3 = 7*(i) + 7003;
      A3_1_3 = 7*(i) + 7004;
      A3_2_3 = 7*(i) + 7005;
      A3_3_3 = 7*(i) + 7006;
      A4_3 = 7*(i) + 7007;
    
      // Fourth batch
    
      A7_4 = 7*(i)  + 10501;
      A6_1_4 = 7*(i) + 10502;
      A6_2_4 = 7*(i) + 10503;
      A3_1_4 = 7*(i) + 10504;
      A3_2_4 = 7*(i) + 10505;
      A3_3_4 = 7*(i) + 10506;
      A4_4 = 7*(i) + 10507;
    
      //  --- C12E2-M --- //
    
    
      A13 = 7*(i) + 1751;
      A12_1= 7*(i)+ 1752;
      A12_2= 7*(i)+ 1753;
      A9_1 = 7*(i)+ 1754;
      A9_2 = 7*(i)+ 1755;
      A9_3 = 7*(i)+ 1756;
      A10 = 7*(i)+ 1757;
    
      /* Mimic atoms - Second batch */
      
      A13_2 = 7*(i) + 5251;
      A12_1_2 = 7*(i)+ 5252;
      A12_2_2 = 7*(i)+ 5253;
      A9_1_2 = 7*(i)+ 5254;
      A9_2_2 = 7*(i)+ 5255;
      A9_3_2 = 7*(i)+ 5256;
      A10_2 = 7*(i)+ 5257;
    
// Mimic atoms - Third batch 
    
      A13_3 = 7*(i) + 8751;
      A12_1_3 = 7*(i)+ 8752;
      A12_2_3 = 7*(i)+ 8753;
      A9_1_3 = 7*(i)+ 8754;
      A9_2_3 = 7*(i)+ 8755;
      A9_3_3 = 7*(i)+ 8756;
      A10_3 = 7*(i)+ 8757;
    
// Mimic atoms - Fourth batch */
  
      A13_4 = 7*(i) + 12251;
      A12_1_4 = 7*(i)+ 12252;
      A12_2_4 = 7*(i)+ 12253;
      A9_1_4 = 7*(i)+ 12254;
      A9_2_4 = 7*(i)+ 12255;
      A9_3_4 = 7*(i)+ 12256;
      A10_4 = 7*(i)+ 12257;

      
        
      // First batch
      C12E2_struct[i].index[0] = A7;
      C12E2_struct[i].index[1] = A6_1;
      C12E2_struct[i].index[2] = A6_2;
      C12E2_struct[i].index[3] = A3_1;
      C12E2_struct[i].index[4] = A3_2;
      C12E2_struct[i].index[5] = A3_3;
      C12E2_struct[i].index[6] = A4;

      // Second batch
      C12E2_struct[i+250].index[0] = A7_2;
      C12E2_struct[i+250].index[1] = A6_1_2;
      C12E2_struct[i+250].index[2] = A6_2_2;
      C12E2_struct[i+250].index[3] = A3_1_2;
      C12E2_struct[i+250].index[4] = A3_2_2;
      C12E2_struct[i+250].index[5] = A3_3_2;
      C12E2_struct[i+250].index[6] = A4_2;
      
      // Third batch 
      C12E2_struct[i+500].index[0] = A7_3;
      C12E2_struct[i+500].index[1] = A6_1_3;
      C12E2_struct[i+500].index[2] = A6_2_3;
      C12E2_struct[i+500].index[3] = A3_1_3;
      C12E2_struct[i+500].index[4] = A3_2_3;
      C12E2_struct[i+500].index[5] = A3_3_3;
      C12E2_struct[i+500].index[6] = A4_3;
      
      // Fourth batch
      C12E2_struct[i+750].index[0] = A7_4;
      C12E2_struct[i+750].index[1] = A6_1_4;
      C12E2_struct[i+750].index[2] = A6_2_4;
      C12E2_struct[i+750].index[3] = A3_1_4;
      C12E2_struct[i+750].index[4] = A3_2_4;
      C12E2_struct[i+750].index[5] = A3_3_4;
      C12E2_struct[i+750].index[6] = A4_4;
      
      // Mimic assignment
      
      
      // First batch
      C12E2M_struct[i].index[0] = A13;
      C12E2M_struct[i].index[1] = A12_1;
      C12E2M_struct[i].index[2] = A12_2;
      C12E2M_struct[i].index[3] = A9_1;
      C12E2M_struct[i].index[4] = A9_2;
      C12E2M_struct[i].index[5] = A9_3;
      C12E2M_struct[i].index[6] = A10;
      
      // Second batch
      C12E2M_struct[i+250].index[0] = A13_2;
      C12E2M_struct[i+250].index[1] = A12_1_2;
      C12E2M_struct[i+250].index[2] = A12_2_2;
      C12E2M_struct[i+250].index[3] = A9_1_2;
      C12E2M_struct[i+250].index[4] = A9_2_2;
      C12E2M_struct[i+250].index[5] = A9_3_2;
      C12E2M_struct[i+250].index[6] = A10_2;
      
      // Third batch 
      C12E2M_struct[i+500].index[0] = A13_3;
      C12E2M_struct[i+500].index[1] = A12_1_3;
      C12E2M_struct[i+500].index[2] = A12_2_3;
      C12E2M_struct[i+500].index[3] = A9_1_3;
      C12E2M_struct[i+500].index[4] = A9_2_3;
      C12E2M_struct[i+500].index[5] = A9_3_3;
      C12E2M_struct[i+500].index[6] = A10_3;
      
      // Fourth batch
      C12E2M_struct[i+750].index[0] = A13_4;
      C12E2M_struct[i+750].index[1] = A12_1_4;
      C12E2M_struct[i+750].index[2] = A12_2_4;
      C12E2M_struct[i+750].index[3] = A9_1_4;
      C12E2M_struct[i+750].index[4] = A9_2_4;
      C12E2M_struct[i+750].index[5] = A9_3_4;
      C12E2M_struct[i+750].index[6] = A10_4; 	    
    }
  }
*/
