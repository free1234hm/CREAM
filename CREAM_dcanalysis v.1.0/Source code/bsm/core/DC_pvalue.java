package bsm.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import org.apache.commons.math3.distribution.NormalDistribution;
import bsm.BSM_DataSet;
import ssclust.Pearsonr;

public class DC_pvalue {
	
	public double[] gcpathway (BSM_DataSet controldata, BSM_DataSet casedata, List<Integer> list, int samplenum){
		Pearsonr r = new Pearsonr();
		double[] result = new double[3];
		double[] listcontrol = new double[list.size()*(list.size()-1)/2];
		double[] listcase = new double[list.size()*(list.size()-1)/2];
		double controlsum = 0;
		double casesum = 0;
		int permutation1 = 30;
		
		Random rd = new Random();
		if(controldata.numcols > samplenum) {
			for(int ran=0;ran<permutation1;ran++) {
				List<Integer> controlsample = new ArrayList<Integer>();
				while(controlsample.size()<samplenum) {
					int num = rd.nextInt(controldata.numcols);
					if(!controlsample.contains(num)) {
						controlsample.add(num);
					}
				}
				int count = 0;
				for(int i=0;i<list.size();i++){
					double[] control1 = new double[samplenum];
					for(int j=0;j<samplenum;j++) {
						control1[j] = controldata.controlnorm[list.get(i)][controlsample.get(j)];
					}
					for(int j=i+1;j<list.size();j++) {
						double[] control2 = new double[samplenum];
						for(int m=0;m<samplenum;m++) {
							control2[m] = controldata.controlnorm[list.get(j)][controlsample.get(m)];
						}
						double controlR = Math.abs(r.pearsonr(control1, control2));
						listcontrol[count] += controlR;
						controlsum += controlR;
						count++;
					}
				}
			}
			controlsum = controlsum/listcontrol.length/permutation1;
			for(int i=0;i<listcontrol.length;i++) {
				listcontrol[i] = listcontrol[i]/permutation1;
			}
		} else {
			int count = 0;
			for(int i=0;i<list.size();i++){
				double[] control1 = controldata.controlnorm[list.get(i)];
				for(int j=i+1;j<list.size();j++) {
					double[] control2 = controldata.controlnorm[list.get(j)];
					double controlR = Math.abs(r.pearsonr(control1, control2));
					listcontrol[count] = controlR;
					controlsum += controlR;
					count++;
				}
			}
			controlsum = controlsum/count;
		}
		
		if(casedata.numcols > samplenum) {
			for(int ran=0;ran<permutation1;ran++) {
				List<Integer> casesample = new ArrayList<Integer>();
				while(casesample.size()<samplenum) {
					int num = rd.nextInt(casedata.numcols);
					if(!casesample.contains(num)) {
						casesample.add(num);
					}
				}
				int count = 0;
				for(int i=0;i<list.size();i++){
					double[] case1 = new double[samplenum];
					for(int j=0;j<samplenum;j++) {
						case1[j] = casedata.controlnorm[list.get(i)][casesample.get(j)];
					}
					for(int j=i+1;j<list.size();j++) {
						double[] case2 = new double[samplenum];
						for(int m=0;m<samplenum;m++) {
							case2[m] = casedata.controlnorm[list.get(j)][casesample.get(m)];
						}
						double caseR = Math.abs(r.pearsonr(case1, case2));
						listcase[count] += caseR;
						casesum += caseR;
						count++;
					}
				}
			}
			casesum = casesum/listcase.length/permutation1;
			for(int i=0;i<listcase.length;i++) {
				listcase[i] = listcase[i]/permutation1;
			}
		} else {
			int count = 0;
			for(int i=0;i<list.size();i++){
				double[] case1 = casedata.controlnorm[list.get(i)];
				for(int j=i+1;j<list.size();j++) {
					double[] case2 = casedata.controlnorm[list.get(j)];
					double caseR = Math.abs(r.pearsonr(case1, case2));
					listcase[count] = caseR;
					casesum += caseR;
					count++;
				}
			}
			casesum = casesum/count;
		}
		
		result[0] = controlsum;
		result[1] = casesum;
		List<Double> listall = new ArrayList<Double>();
		
			double DS=0;
			int num = listcontrol.length;
			for(int i=0;i<listcontrol.length;i++) {
				listall.add(listcontrol[i]);
			}
			for(int i=0;i<listcase.length;i++) {
				listall.add(listcase[i]);
			}
			
			for(int m=0;m<num;m++){
				DS = DS + (listall.get(m)-listall.get(m+num)) / num;
			}
			//DS=Math.sqrt(DS);
			int permutation = 1000;
			double Sum = 0;
			double Sumsq = 0;
			for(int random=0;random<permutation;random++){
				double ranDS=0;
				Collections.shuffle(listall);
				for(int m=0;m<num;m++){
					ranDS = ranDS + (listall.get(m)-listall.get(m+num)) / num;
				}
				//ranDS=Math.sqrt(ranDS);
				Sum += ranDS;
				Sumsq += Math.pow(ranDS, 2);
			}
			double mean = Sum / permutation;
			double sigma = Math.sqrt(Sumsq/permutation - Math.pow(mean,2));
			NormalDistribution nor = new NormalDistribution(mean, sigma);
			result[2] = 1-nor.cumulativeProbability(DS);
			System.out.println(mean+"\t"+sigma+"\t"+DS+"\t"+result[2]);
		return result;
	}
	
}
