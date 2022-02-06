time = [50:600]; % [s]
temperature = 950 - 932*(0.996716).^time;

plot(time,temperature,'LineWidth',1.5)
xticks([50:100:550])
xlim([50,600])
xticklabels([60,50,40,30,20,10])
xlabel("Velocity of Plates [mm/s]")
ylabel("Temperature at Exit [C]")
title("Temperature of Plates at Exit of Oven")